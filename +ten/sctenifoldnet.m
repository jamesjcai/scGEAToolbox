function [T, A0, A1] = sctenifoldnet(X0, X1, genelist, varargin)
% T=sctenifoldnet_m(X0,X1,genelist);
%
% X0 and X1 are gene x cell matrices
%
% Name-value options: 'qqplot', 'smplmethod', 'tdmethod', 'nsubsmpl',
% 'csubsmpl', 'savegrn', 'useparallel'.
%
% 'useparallel' (default false) runs TEN.I_NC's network-construction loop as a
% parfor, one worker per subsampled network. See I_NC's header before turning
% it on: the serial path is already parallel through multithreaded BLAS, so
% the realistic gain is modest rather than nsubsmpl-fold, and the parallel
% branch draws different bootstrap subsamples than the serial one - equivalent
% in distribution, reproducible from a seed, but not the same cells.
%
import ten.*

if ~(ismcc || isdeployed)
    if exist(['@tensor', filesep, 'tensor.m'], 'file') ~= 2
        if ispref('scgeatoolbox', 'tensor_toolbox_path')
            pth = getpref('scgeatoolbox', 'tensor_toolbox_path');
            if isfolder(pth)
                addpath(pth);
            else
                error('sctenifoldnet:missingToolbox', ...
                    'Tensor Toolbox path not found: %s\nRe-install via Tools > Install Tensor Toolbox.', pth);
            end
        else
            error('sctenifoldnet:missingToolbox', ...
                'Tensor Toolbox is not installed. Install it via Tools > Install Tensor Toolbox.');
        end
    end
end

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
parse(p, varargin{:});
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

validg = ~ismember(upper(genelist), upper(pkg.i_get_ribosomalgenes));
genelist = genelist(validg);
X0 = sc_norm(X0(validg, :), "type", "libsize");
X1 = sc_norm(X1(validg, :), "type", "libsize");

X0 = log1p(X0);
X1 = log1p(X1);

% X0=sc_transform(X0);
% X1=sc_transform(X1);

rng('default');

tic
disp('Sample 1/2 ...')
[XM] = i_nc(X0, nsubsmpl, 3, csubsmpl, usebootstrp, useparallel);
toc
tic
disp('Tensor decomposition')
[A0] = i_td1(XM, tdmethod);
% [A0]=e_filtadjc(A0,0.95);
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
% A1=e_filtadjc(A1,0.95);
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
