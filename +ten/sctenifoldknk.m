function [T, A0, genelist] = sctenifoldknk(X, genelist, kogene, varargin)
% T=sctenifoldknk(X,genelist,"Foxp3");
%
% X is a gene x cell matrix from the wild-type
%
% The knockout removes KOGENE's outgoing edges, as in R; see TEN.I_KNK for the
% 2026-09-25 fix of its direction. A0 is returned in NET.PCRNET's orientation
% (row k = regression of gene k), over the genes of the third output.
%
% THE TOOLBOX VARIANT IS THE DEFAULT ('reference', false): libsize + log1p,
% no QC, Tensor Toolbox CP of rank NSUBSMPL, normalized Laplacian, and at
% least 10 outgoing links on the KO gene. The third output is then GENELIST.
%
% WHY THE TOOLBOX VARIANT IS THE DEFAULT. Reference mode was the default for
% part of 2026-09-25 and was switched back after a paired benchmark
% (pipelines/tenifold_laplacian, README "Why reference mode scores lower"):
% on planted co-expression changes R's pipeline scored AUROC 0.72 against
% 0.93 for the toolbox variant, the loss coming from R's rounding of the denoised network to one
% decimal, its rank-3 decomposition and CPM without log1p. Use reference mode
% when the numbers have to agree with R or scTenifoldpy.
%
% REFERENCE MODE ('reference', true) reproduces scTenifoldKnk() of the R
% package (1.1.4) run with its defaults as scTenifoldKnk(X, gKO, seed = SEED):
%
%     [T, A0, g] = ten.sctenifoldknk(X, genelist, "Foxp3", 'reference', true);
%
% X must then be raw counts. Every step follows R: TEN.I_SCQC for QC and CPM
% (no log1p; KOGENE must pass QC), TEN.I_NCREF with q = 0.9, TEN.I_CPALS
% with K = 3 rounded to 3 decimals, and TEN.I_KNK in its reference mode
% (unnormalized, unsymmetrized alignment in 2 dimensions; TEN.I_DR reference
% output, sorted by p). 'seed' (default 1) is R's seed; 'smplmethod' and
% 'tdmethod' are ignored, and so is 'sorttable': the table is always sorted,
% as R's is. A0 has a zero diagonal. G, the genes after QC, is what A0 and T
% are over - fewer than GENELIST when QC drops any, so index A0 with G, not
% with GENELIST. Results match R to floating-point precision, except for
% edges tied at the quantile threshold or values right on a rounding
% boundary.
%
% CHECK THE RESULT AGAINST A NULL BEFORE USING IT. T's p-values come from
% TEN.I_DR, which fits a chi-square to drdist/mean(drdist). That asks whether
% a gene moved more than the average gene moved in this run - not whether it
% moved because of KOGENE. On a network with strong hub structure the same
% genes top the ranking whichever row is zeroed, and the enrichment computed
% from them describes the network rather than the knockout.
%
% TEN.KNKNULLCONTROL settles it by knocking out random genes from the same A0
% and rescoring against that background. It costs minutes, needs the A0
% returned here, and reports a verdict:
%
%     [T, A0, g] = ten.sctenifoldknk(X, genelist, "Foxp3");
%     [Tnull, S] = ten.knknullcontrol(A0, "Foxp3", g);
%     disp(S.Verdict)
%
% see also: TEN.KNKNULLCONTROL, TEN.I_KNK, TEN.I_DR, TEN.SCTENIFOLDNET
import ten.*

if nargin < 3
    error(sprintf('USAGE: T=sctenifoldknk(X,genelist,kogene);\n       T=sctenifoldnet_m(X0,X1,genelist,''qqplot'',true);'));
end
if isscalar(kogene) && isnumeric(kogene)
    idx = kogene;
else
    idx = find(genelist == kogene, 1);
    if isempty(idx)
        error("KOGENE should be a member of GENELIST.");
    end
end
p = inputParser;
addOptional(p, 'sorttable', false, @islogical);
addOptional(p, 'smplmethod', "bootstrap", @(x) (isstring(x) | ischar(x)) & ismember(lower(string(x)), ["jackknife", "bootstrap"]));
addOptional(p, 'tdmethod', "CP", @(x) (isstring(x) | ischar(x)) & ismember(upper(string(x)), ["CP", "TUCKER"]));
addOptional(p, 'nsubsmpl', 10, @(x) fix(x) == x & x > 0);
addOptional(p, 'csubsmpl', 500, @(x) fix(x) == x & x > 0);
addOptional(p, 'savegrn', false, @islogical);
addOptional(p, 'reference', false, @islogical);
addOptional(p, 'seed', 1, @(x) isnumeric(x) && isscalar(x) && fix(x) == x);
parse(p, varargin{:});
if p.Results.reference
    [T, A0, genelist] = i_reference(X, genelist, idx, p.Results);
    return
end
dosort = p.Results.sorttable;
tdmethod = p.Results.tdmethod;
nsubsmpl = p.Results.nsubsmpl;
csubsmpl = p.Results.csubsmpl;
smplmethod = p.Results.smplmethod;
savegrn = p.Results.savegrn;

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

if size(X, 1) ~= length(genelist)
    error('Length of genelist should be the same as the number of rows of X0 or X1.');
end

% Add the Tensor Toolbox from the saved preference if it is not already on
% the path, then error only if it genuinely is not installed.
%
% This used to error outright, which made the function fail in any MATLAB
% that had not happened to call ten.sctenifoldnet first - that one does the
% addpath itself, so an interactive session picked the path up as a side
% effect and a fresh `matlab -batch` did not. A long unattended run would
% then die at its first contrast. ten.check_tensor_toolbox is the toolbox's
% own helper for exactly this and is what ten.sctenifoldnet's own preamble
% amounts to.
ten.check_tensor_toolbox;

X = sc_norm(X, "type", "libsize");
X = log1p(X);

[XM] = i_nc(X, nsubsmpl, 3, csubsmpl, usebootstrp);
[A0] = i_td1(XM, tdmethod);
    if savegrn
        tstr = matlab.lang.makeValidName(strrep(sprintf("GRN_created_on_% s", datetime)," ", "_at_"));
        save(tstr, 'A0', 'genelist', '-v7.3');
        fprintf('\nConstructed gene regulatory network (GRN) is saved in %s.mat\n', tstr);
    end
T = ten.i_knk(A0, idx, genelist, dosort);
end


function [T, A0, genelist] = i_reference(X, genelist, idx, opts)
% scTenifoldKnk() of the R package; see REFERENCE MODE in the header
if size(X, 1) ~= numel(genelist)
    error('Length of genelist should be the same as the number of rows of X.');
end
genelist = string(genelist(:));
kogene = genelist(idx);
[X, genelist] = ten.i_scqc(X, genelist);
idx = find(genelist == kogene, 1);
if isempty(idx)
    error('ten:sctenifoldknk:koGeneRemoved', ...
        'The KO gene %s is not present in the count matrix after quality control.', kogene);
end

XM = ten.i_ncref(X, 'NumNets', opts.nsubsmpl, 'NumCells', opts.csubsmpl, ...
    'NumComp', 3, 'Q', 0.9, 'Seed', opts.seed);
A0 = ten.i_cpals(XM, 3, NumDecimals=3, Seed=opts.seed);
clear XM
A0(1:(size(A0, 1) + 1):end) = 0;
if opts.savegrn
    tstr = matlab.lang.makeValidName(string(datetime));
    save(sprintf('A0_%s', tstr), 'A0', 'genelist', '-v7.3');
end
T = ten.i_knk(A0, idx, genelist, true, 0, true);
end
