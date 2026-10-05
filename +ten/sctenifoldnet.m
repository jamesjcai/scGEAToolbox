function [T, A0, A1, genelist] = sctenifoldnet(X0, X1, genelist, varargin)
% T=sctenifoldnet_m(X0,X1,genelist);
%
% X0 and X1 are gene x cell matrices
%
% The fourth output is the gene list the returned networks are over. It is
% NOT the list that went in: the ribosomal filter (the default toolbox
% variant) or QC (reference mode) drops genes, so A0 and A1 are
% smaller than X0 and every row after the first dropped gene is shifted.
% Anything that indexes A0 or A1 by position has to use this list --
% T.genelist will not do, because I_DR sorts T.
%
% Name-value options: 'reference', 'seed', 'qqplot', 'nsubsmpl', 'csubsmpl',
% 'savegrn'; and, for the toolbox variant only, 'smplmethod', 'tdmethod',
% 'useparallel'.
%
% THE TOOLBOX VARIANT IS THE DEFAULT ('reference', false): ribosomal genes
% dropped, libsize + log1p, no QC, MATLAB's RNG, TEN.I_NC networks, Tensor
% Toolbox CP of rank NSUBSMPL, normalized Laplacian. The options and notes
% from here to REFERENCE MODE describe it. 'reference', true reproduces the R
% package instead; see REFERENCE MODE below.
%
% WHY THE TOOLBOX VARIANT IS THE DEFAULT. Reference mode was the default for
% part of 2026-09-25 and was switched back after a paired benchmark
% (pipelines/tenifold_laplacian, README "Why reference mode scores lower"):
% on planted co-expression changes R's pipeline scored AUROC 0.72 against
% 0.93 here, the loss coming from R's rounding of the denoised network to one
% decimal, its rank-3 decomposition and CPM without log1p. Use reference mode
% when the numbers have to agree with R or scTenifoldpy.
%
% 'useparallel' (default false) runs TEN.I_NC's network-construction loop as a
% parfor, one worker per subsampled network. See I_NC's header before turning
% it on: the serial path is already parallel through multithreaded BLAS, so
% the realistic gain is modest rather than nsubsmpl-fold, and the parallel
% branch draws different bootstrap subsamples than the serial one - equivalent
% in distribution, reproducible from a seed, but not the same cells.
%
% Manifold alignment uses the normalized Laplacian since 2026-09-24; see
% TEN.I_MA for why, and for how to reproduce earlier results.
%
% REFERENCE MODE ('reference', true) reproduces scTenifoldNet() of the R
% package (1.4.3) run with its defaults as scTenifoldNet(X0, X1, seed = SEED):
%
%     [T, A0, A1, g] = ten.sctenifoldnet(X0, X1, genelist, 'reference', true);
%
% X0 and X1 must be raw counts. Every step follows R instead of the toolbox
% variant:
%   QC + CPM     TEN.I_SCQC (scQC, cpmNormalization); no log1p, ribosomal
%                genes kept; the genes are those passing QC in both samples
%   networks     TEN.I_NCREF: R's sample() through TEN.RRANDOM, q = 0.95 on
%                R's quantile, max scaling; 'nsubsmpl' and 'csubsmpl' apply
%   tensor       TEN.I_CPALS: R's CP-ALS, K = 3, rounded to 1 decimal; no
%                Tensor Toolbox needed
%   alignment    unnormalized Laplacian, 30 dimensions
%   regulation   TEN.I_DR reference mode: adds distance and Box-Cox Z,
%                noise-level distances get p = 1, sorted by p
% 'seed' (default 1) seeds the subsampling and the decomposition as R's
% seed does; 'smplmethod', 'tdmethod' and 'useparallel' are ignored. A0 and
% A1 are the tensor networks with a zero diagonal, as R returns them. The
% results match R to floating-point precision (the manifold up to the sign
% of each column); an edge whose weight ties the quantile threshold, or a
% network value right on a rounding boundary, can still come out
% differently.
%
import ten.*

if nargin < 2
    error(sprintf('USAGE: T=sctenifoldnet(X0,X1);\n       T=sctenifoldnet(X0,X1,genelist,''qqplot'',true);'));
end
if nargin < 3, genelist = string(num2cell(1:size(X0, 1)))'; end

p = inputParser;
addOptional(p, 'qqplot', false, @islogical);
addOptional(p, 'smplmethod', "Bootstrap", @(x) (isstring(x) | ischar(x)) & ismember(lower(string(x)), ["jackknife", "bootstrap"]));
addOptional(p, 'tdmethod', "CP", @(x) (isstring(x) | ischar(x)) & ismember(upper(string(x)), ["CP", "TUCKER"]));
addOptional(p, 'nsubsmpl', 10, @(x) fix(x) == x & x > 0);
addOptional(p, 'csubsmpl', 500, @(x) fix(x) == x & x > 0);
addOptional(p, 'savegrn', true, @islogical);
addOptional(p, 'useparallel', false, @islogical);
addOptional(p, 'reference', false, @islogical);
addOptional(p, 'seed', 1, @(x) isnumeric(x) && isscalar(x) && fix(x) == x);
parse(p, varargin{:});

if p.Results.reference
    [T, A0, A1, genelist] = i_reference(X0, X1, genelist, p.Results);
    if p.Results.qqplot
        figure;
        e_mkqqplot(T);
    end
    return
end

if ~(ismcc || isdeployed)
    if exist(['@tensor', filesep, 'tensor.m'], 'file') ~= 2
        if ispref('scgeatoolbox', 'tensor_toolbox_path')
            pth = getpref('scgeatoolbox', 'tensor_toolbox_path');
            if isfolder(pth)
                addpath(pth);
            else
                error('sctenifoldnet:missingToolbox', ...
                    'Tensor Toolbox path not found: %s\nRe-install via Setup > Install Tensor Toolbox.', pth);
            end
        else
            error('sctenifoldnet:missingToolbox', ...
                'Tensor Toolbox is not installed. Install it via Setup > Install Tensor Toolbox.');
        end
    end
end

doqqplot = p.Results.qqplot;
tdmethod = p.Results.tdmethod;
nsubsmpl = p.Results.nsubsmpl;
csubsmpl = p.Results.csubsmpl;
smplmethod = p.Results.smplmethod;
savegrn = p.Results.savegrn;
useparallel = p.Results.useparallel;

switch upper(tdmethod)
    case "CP"
        tdmethod = 1;
    case "TUCKER"
        tdmethod = 2;
end
switch lower(smplmethod)
    case "jackknife"
        usebootstrp = false;
    case "bootstrap"
        usebootstrp = true;
end

if size(X0, 1) ~= size(X1, 1)
    error('X0 and X1 need the same number of rows.');
end
if size(X0, 1) ~= length(genelist)
    error('Length of genelist should be the same as the number of rows of X0 or X1.');
end

if exist('tensor.m', 'file') ~= 2
    error('Need tensor_toolbox');
end
if isempty(which('net.pcrnet'))
    error('Need net.pcrnet in scGEAToolbox https://github.com/jamesjcai/scGEAToolbox');
end

% GENELIST is narrowed here, which is why it is also an output: A0 and A1
% below are over the surviving genes only.
validg = ~ismember(upper(genelist), upper(pkg.i_get_ribosomalgenes));
genelist = genelist(validg);
X0 = sc_norm(X0(validg, :), "type", "libsize");
X1 = sc_norm(X1(validg, :), "type", "libsize");

X0 = log1p(X0);
X1 = log1p(X1);

% Save and restore the caller's random stream. Seeding the bootstrap subsamples is
% fine; leaving the session parked on that seed is not -- it then
% governs every later tsne, umap and clustering call in the session.
% ten.sctenifoldnetstability's header documents this determinism and
% relies on it, so the reset itself stays.
rngState = rng();
restoreRng = onCleanup(@() rng(rngState));
rng('default');

tic
disp('Sample 1/2 ...')
[XM] = i_nc(X0, nsubsmpl, 3, csubsmpl, usebootstrp, useparallel);
toc
tic
disp('Tensor decomposition')
[A0] = i_td1(XM, tdmethod);
toc
if savegrn
    tic
    disp('Saving scGRN network 1/2')
    tstr = matlab.lang.makeValidName(string(datetime));
    save(sprintf('A0_%s', tstr), 'A0', 'genelist', '-v7.3');
    toc
end
tic
disp('Sample 2/2 ...')
[XM] = i_nc(X1, nsubsmpl, 3, csubsmpl, usebootstrp, useparallel);
toc
tic
disp('Tensor decomposition')
A1 = i_td1(XM, tdmethod);
toc
if savegrn
    tic
    disp('Saving scGRN network 2/2')
    save(sprintf('A1_%s', tstr), 'A1', 'genelist', '-v7.3');
    toc
end

A0sym = 0.5 * (A0 + A0.');
A1sym = 0.5 * (A1 + A1.');
tic
disp('Manifold alignment')
[aln0, aln1] = ten.i_ma(A0sym, A1sym);
toc

T = ten.i_dr(aln0, aln1, genelist);

if doqqplot
    figure;
    e_mkqqplot(T);
end
end


function [T, A0, A1, genelist] = i_reference(X0, X1, genelist, opts)
% scTenifoldNet() of the R package; see REFERENCE MODE in the header
if size(X0, 1) ~= numel(genelist) || size(X1, 1) ~= numel(genelist)
    error('Length of genelist should be the same as the number of rows of X0 and X1.');
end
genelist = string(genelist(:));
[X0, g0] = ten.i_scqc(X0, genelist);
[X1, g1] = ten.i_scqc(X1, genelist);

% intersect(rownames(X), rownames(Y)), in X's order
genelist = unique(g0(ismember(g0, g1)), 'stable');
[~, i0] = ismember(genelist, g0);
[~, i1] = ismember(genelist, g1);
X0 = X0(i0, :);
X1 = X1(i1, :);
fprintf('Shared genes: %d\n', numel(genelist));

ncArgs = {'NumNets', opts.nsubsmpl, 'NumCells', opts.csubsmpl, ...
    'NumComp', 3, 'Q', 0.95, 'Seed', opts.seed};
disp('Sample 1/2 ...')
XM = ten.i_ncref(X0, ncArgs{:});
disp('Tensor decomposition')
A0 = ten.i_cpals(XM, 3, NumDecimals=1, Seed=opts.seed);
disp('Sample 2/2 ...')
XM = ten.i_ncref(X1, ncArgs{:});
disp('Tensor decomposition')
A1 = ten.i_cpals(XM, 3, NumDecimals=1, Seed=opts.seed);
clear XM

disp('Manifold alignment')
[aln0, aln1] = ten.i_ma(0.5*(A0 + A0.'), 0.5*(A1 + A1.'), 30, "unnormalized");
T = ten.i_dr(aln0, aln1, genelist, true, true);

% The alignment ignores self-loops; R drops them from the returned networks
A0(1:(size(A0, 1) + 1):end) = 0;
A1(1:(size(A1, 1) + 1):end) = 0;
if opts.savegrn
    tstr = matlab.lang.makeValidName(string(datetime));
    save(sprintf('A0_%s', tstr), 'A0', 'genelist', '-v7.3');
    save(sprintf('A1_%s', tstr), 'A1', 'genelist', '-v7.3');
end
end
