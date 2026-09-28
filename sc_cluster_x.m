function [c] = sc_cluster_x(X, k, varargin)
% sc_cluster_x - cluster cells using UMI matrix X
%
% see also: sc_cluster_s

p = inputParser;
defaultType = 'sc3';
% 'specter' was removed: run.ml_Specter split 160 cells 12/148 where the
% truth was 80/80 -- ARI 0.014 on data plain k-means clusters perfectly,
% identically across six seeds. It is retired to unused/.
validTypes = {'sc3', 'simlr', 'soptsc', 'sinnlrr'};
checkType = @(x) any(validatestring(x, validTypes));

checkK = @(x) (x > 0) && isnumeric(x) && isscalar(x);

addRequired(p, 'X', @isnumeric);
addRequired(p, 'k', checkK);
addOptional(p, 'type', defaultType, checkType);
addOptional(p, 'usehvgs', true);

parse(p, X, k, varargin{:});

if p.Results.usehvgs
    disp('Using 2000 HVGs.')
    % SC_ANALYTICFIT scores each gene against the closed-form
    % gamma-Poisson curve implied by the library sizes, the same
    % ranker SINGLECELLEXPERIMENT.EMBEDCELLS uses. It returns gene
    % names rather than a sorted matrix, so the rows are picked out
    % by name. The names are synthetic and unique here, and gene
    % order does not matter to any of the clustering backends.
    g = "gene_" + string((1:size(X, 1)).');
    T = sc_analyticfit(X, g);
    nkeep = min(height(T), 2000);
    X = X(ismember(g, T.genes(1:nkeep)), :);
end

switch p.Results.type
    case 'simlr'
        [c] = run.ml_SIMLR(X, k, true);
    case 'soptsc'
        % Symmetric NMF for cell clustering
        % https://www.biorxiv.org/content/biorxiv/early/2019/01/01/168922.full.pdf
        % disp('To specify k, use RUN_SOPTSC(X,''k'',k).');
        [c] = run.ml_SoptSC(X, 'k', k, 'donorm', true);
    case 'sc3'
        [c] = run.ml_SC3(X, k);
    case 'sinnlrr'
        [c] = run.ml_SinNLRR(X, k);
end
end