function X = ml_MAGIC(X, donorm)
%ML_MAGIC  MAGIC imputation (van Dijk et al. 2018), genes-by-cells in and out.
%
%   X = RUN.ML_MAGIC(X) library-size normalises X and returns the imputed
%   counts. RUN.ML_MAGIC(X, false) skips the normalisation.
%
%   Needs at least MINCELLS cells; see the check below for why that is an
%   error rather than something to work around.

if nargin < 2, donorm = true; end

% SIZE CHECK, in the caller's genes-by-cells orientation.
%
% RANDPCA -- which RUN_MAGIC_ORIGINAL calls as randPCA(data', npca) --
% refuses more components than the smallest dimension of its input, and
% this wrapper passed a hard-coded npca of 100. So any matrix with fewer
% than 100 genes OR fewer than 100 cells failed, and failed opaquely:
% sc_impute on a 60-by-300 matrix came back with "Input 2 must be <= the
% smallest dimension of Input 1", which names neither npca nor the
% caller's data. Measured: 150x300 worked; 60x300, 150x60 and 40x40 all
% failed that way.
%
% The two dimensions are not symmetric, so they are handled differently.
% Genes only bound npca, so npca is clamped below and a gene-poor matrix
% now works. Cells drive the k = 15 neighbour graph and the optimal
% diffusion-time search, and too few of them does not error -- it fails to
% terminate. So cells are checked here, up front, and never clamped into
% that. Erroring immediately with a legible message is the point; the
% opaque failure it replaces was at least prompt.
[nGenesIn, nCellsIn] = size(X);
minCells = 100;
if nCellsIn < minCells
    error('run:ml_MAGIC:tooFewCells', ...
        ['MAGIC needs at least %d cells to build its k-nearest-neighbour ', ...
        'graph and search for a diffusion time; this matrix has %d. Use ', ...
        'more cells, or a method that does not diffuse over a cell graph.'], ...
        minCells, nCellsIn);
end
if nGenesIn < 2
    error('run:ml_MAGIC:tooFewGenes', ...
        'MAGIC needs at least two genes; this matrix has %d.', nGenesIn);
end

pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '..', 'external', 'ml_MAGIC');
if ~(ismcc || isdeployed)
    addpath(pth);
end
% gene_names=cellstr(gl123);  MAGIC needs [cells x genes]
% data=X';

% library size normalization
% libsize = sum(data,2);
% data = bsxfun(@rdivide, data, libsize) * median(libsize);
if donorm
    X = sc_norm(X, 'type', 'libsize');
    % log transform -- usually one would log transform the data. Here we don't do it.
    % data = log(data + 0.1);
end

% MAGIC needs [cells x genes]
X = X';

% X is cells-by-genes here, and randPCA sees its transpose, so the bound is
% min(nGenes, nCells) either way. The cell count is already known to be at
% least minCells, so this only ever clamps on gene-poor input.
% The bundled external/ml_MAGIC does a bare "warning off" with
% nothing to undo it (external/ml_MAGIC/randPCA.m line 103), so one call left warnings
% disabled for the rest of the MATLAB session. That is not cosmetic: it
% silences every later warning the user relies on, and it made eleven
% unrelated verifyWarning tests in this repository's suite fail, nowhere
% near the culprit. The same defect in external/ml_SinNLRR was fixed the
% same way in 0827873; +pkg/e_randPCA.m, the in-toolbox copy of this very
% file, has the offending line commented out. Restoring here rather than
% editing the third-party file keeps the fix in code we own.
warnState = warning();
restoreWarn = onCleanup(@() warning(warnState));
npca = min(100, min(size(X)));

[pc_imputed, U, ~] = run_magic_original(X, 'npca', npca, 'k', 15, 'a', 15, 'make_plot_opt_t', false);

% plot_genes = {'Cdh1', 'Vim', 'Fn1', 'Zeb1'};
% [M_imputed, genes_found] = project_genes(plot_genes, gene_names, pc_imputed, U);
M = pc_imputed * U'; % project
X = M';

end
