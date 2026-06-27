function [A, P, Tmarkers] = r_cogaps(X, g, nPatterns, nIterations, sparseOpt, wkdir, isdebug)
% R_COGAPS - Bayesian NMF pattern discovery via CoGAPS (Bioconductor)
% [A, P, Tmarkers] = run.r_cogaps(X, g, nPatterns, nIterations)
%   X           - gene-by-cell expression matrix
%   g           - gene names (string array, length = size(X,1))
%   nPatterns   - number of patterns to learn (default 8)
%   nIterations - number of MCMC iterations (default 1000)
%   sparseOpt   - use sparse optimization (default true)
% Returns:
%   A        - genes-by-patterns feature loadings (A matrix)
%   P        - cells-by-patterns sample factors (P matrix)
%   Tmarkers - PatternMarker table (gene, pattern, rank)
%
% Reference: Sherman et al., Nat Protoc 18:3690-3731 (2023) [PMID:37828301]

if nargin < 3 || isempty(nPatterns), nPatterns = 8; end
if nargin < 4 || isempty(nIterations), nIterations = 1000; end
if nargin < 5 || isempty(sparseOpt), sparseOpt = true; end
if nargin < 6, wkdir = pkg.i_tempdirfile(); end
if nargin < 7, isdebug = false; end

A = [];
P = [];
Tmarkers = [];

oldpth = pwd();
cleanupCwd = onCleanup(@() cd(oldpth));
[isok, msg, codepth] = commoncheck_R('R_CoGAPS');
if ~isok, error(msg); end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'input.mat', 'genes.txt', 'output.h5', 'output_markers.csv'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

if issparse(X), X = full(X); end
nPatterns = double(nPatterns);
nIterations = double(nIterations);
sparseOpt = logical(sparseOpt);
save('input.mat', 'X', 'nPatterns', 'nIterations', 'sparseOpt', '-v7.3');
writematrix(g(:), 'genes.txt');

Rpath = getpref('scgeatoolbox', 'rexecutablepath', []);
if isempty(Rpath)
    error('R environment has not been set up.');
end

codefullpath = fullfile(codepth, 'script.R');
pkg.i_addwd2script(codefullpath, wkdir, 'R');
pkg.i_runrcode(codefullpath, Rpath);

if exist('output.h5', 'file')
    A = h5read('output.h5', '/A');
    P = h5read('output.h5', '/P');
    % Guard against HDF5 row/column-major transposition between R and MATLAB.
    if size(A, 1) ~= numel(g) && size(A, 2) == numel(g)
        A = A';
    end
    if size(P, 1) ~= size(X, 2) && size(P, 2) == size(X, 2)
        P = P';
    end
end
if exist('output_markers.csv', 'file')
    Tmarkers = readtable('output_markers.csv', 'Delimiter', ',');
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
end
