function [X] = sc_transform(X, varargin)
% SC_TRANSFORM Transformations for single-cell data
%
% Supports:
%   - PearsonResiduals
%   - kNNSmoothing
%   - SCTransform / SCT      (R engine, Seurat::SCTransform, vst.flavor "v2")
%   - SCTransformMATLAB      (native MATLAB port, no R required)
%   - FreemanTukey
%
% The SCTransform and SCTransformMATLAB options implement the same
% algorithm; the former calls R/Seurat, the latter is a pure-MATLAB
% reproduction (see sc_sctransformv2) that agrees with R to a correlation
% of ~0.99999 on the Pearson residuals.
%
% Example: X = sc_transform(X, 'type', 'PearsonResiduals');

% https://www.biorxiv.org/content/10.1101/2021.06.24.449781v1.full
% acosh transformation based on the delta method
% shifted logarithm (log(x + c)) with a pseudo-count c, so that it approximates the acosh transformation
% randomized quantile and Pearson residuals

p = inputParser;
defaultType = 'PearsonResiduals';
validTypes = {'PearsonResiduals', 'kNNSmoothing', 'FreemanTukey', ...
    'SCTransform', 'SCT', 'SCTransformR', ...          % R (Seurat) engine
    'SCTransformMATLAB', 'SCTransformM'};              % native MATLAB port
checkType = @(x) any(validatestring(x, validTypes));

addRequired(p, 'X', @isnumeric);
addOptional(p, 'type', defaultType, checkType);
parse(p, X, varargin{:});

ptype = lower(p.Results.type);
% PearsonResiduals reads X block by block as it is, sparse or not; the
% other methods take it dense.
if issparse(X) && ~strcmp(ptype, 'pearsonresiduals'), X = full(X); end

switch ptype
    case 'pearsonresiduals'
        % analytic Pearson residuals
        % https://doi.org/10.1101/2020.12.01.405886
        % https://gist.github.com/hypercompetent/51a3c428745e1c06d826d76c3671797c
        X = i_pearsonresiduals(X);
    case 'knnsmoothing'
        % K-nearest neighbor smoothing for high-throughput single-cell RNA-Seq data
        % https://doi.org/10.1101/217737
        %
        X = knn_smooth(X, 5, 10);

    case {'sct', 'sctransform', 'sctransformr'}
        % sctransform via R/Seurat (SCTransform, vst.flavor "v2")
        % sctransform: Variance Stabilizing Transformations for Single Cell UMI Data
        % Hafemeister & Satija 2019
        % https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
        % https://satijalab.org/seurat/archive/v4.3/sctransform_v2_vignette
        [X] = run.r_SeuratSctransform(X, string(1:size(X, 1)));

    case {'sctransformmatlab', 'sctransformm'}
        % sctransform via the native MATLAB port (no R required).
        % Reproduces sctransform::vst(vst.flavor = "v2") Pearson residuals.
        X = sc_sctransformv2(X);
        % Genes below the min_cells threshold are returned as NaN rows by
        % the reference algorithm; zero them for a usable dense matrix,
        % matching the PearsonResiduals convention above.
        X(isnan(X)) = 0;

    case 'freemantukey'
        % https://github.com/flo-compbio/monet/blob/master/monet/util/expression.py
        % Applies the Freeman-Tukey transformation to stabilize variance."
        % https://www.biorxiv.org/content/10.1101/2020.06.08.140673v2.full
        % https://www.nature.com/articles/nmeth.2930
        X = sc_norm(X, 'type', 'deseq');
        X = sqrt(X) + sqrt(X + 1);
    otherwise
        error('sc_transform:InvalidType', 'Unknown transformation type: %s', p.Results.type);
end
end


function R = i_pearsonresiduals(X)
% Analytic Pearson residuals, (x - u)/sqrt(u + u^2/100) with
% u = rowsum*colsum/total, NaN (where u = 0) set to 0 and clipped to
% +/- sqrt(number of cells).
%
% The result is dense by nature -- a zero count has residual -u/s -- so
% one genes x cells matrix is unavoidable. Everything else is not. This
% used to FULL() the input and then hold U, S and X - U as further dense
% genes x cells arrays: about five at once, over 30 GB at 20000 x 50000.
% Now the output is allocated once and filled a block of cells at a time
% from X as given (a sparse X stays sparse), with the same arithmetic on
% every element.
[G, C] = size(X);
rs = full(sum(X, 2));
cs = full(sum(X, 1));
tot = sum(rs);
sn = sqrt(C);
% U is 0 -- and the residual 0/0 = NaN -- exactly where a gene or a cell
% has no counts at all, so those rows and columns are zeroed directly
% rather than by scanning every element for NaN. MIN/MAX then clip in one
% pass; with no NaN left they give what the two masked assignments did.
zeroRow = rs == 0;
R = zeros(G, C);
step = max(1, floor(1.6e7/max(G, 1)));    % ~128 MB of temporaries per block
for c0 = 1:step:C
    cols = c0:min(C, c0 + step - 1);
    u = (rs * cs(cols)) ./ tot;
    r = (full(X(:, cols)) - u) ./ sqrt(u + (u.^2) ./ 100);
    r(zeroRow, :) = 0;
    r(:, cs(cols) == 0) = 0;
    R(:, cols) = min(max(r, -sn), sn);
end
end


% K-nearest neighbor smoothing for high-throughput scRNA-Seq data
% (Matlab implementation.)

% Author: Maayan Baron <Maayan.Baron@nyumc.org>
% Copyright (c) 2017, 2018 New York University

function [mat_smooth] = knn_smooth(raw_mat, k, varargin)
% This function smooth computes the smoothed matrix using K-nearest neighbor
% by Wagner et al. (2018). The output is a smoothed matrix, not normalized or
% transformed.
% Dependencies: Randomized Singular Value Decomposition (rsvd) function:
% https://www.mathworks.com/matlabcentral/fileexchange/47835-randomized-singular-value-decomposition
%
%                           INPUT ARGUMENTS
%                           ---------------
%   [mat_smooth] = knn_smooth(raw_mat,k) computes the smoothed expression of
%   raw_mat where rows are genes/features and columns are samples/cells using
%   k neighbours.
%   varargin:
%   'num_of_pc' (default = 10) perform PCA before computing distances
%
%                         OUTPUT ARGUMENTS
%                         ----------------
%   The smoothed expression matrix.

% default
num_of_pc = 10;

if ~isempty(varargin)
    num_of_pc = varargin{1};
end

mat_smooth = raw_mat;
num_of_steps = ceil(log2(k + 1));
disp_text = ['number of steps: ', num2str(num_of_steps)];
disp(disp_text)
for s = 1:num_of_steps
    k_step = min(2^s - 1, k);
    mat_tpm = median(sum(mat_smooth)) * bsxfun(@rdivide, mat_smooth, sum(mat_smooth));
    mat_trans = sqrt(mat_tpm) + sqrt(mat_tpm + 1);
    [~, ~] = sort(sum(mat_smooth, 2), 'descend');
    disp_texp = ['preforming pca ', num2str(s), '/', num2str(num_of_steps), ' times'];
    disp(disp_texp)
    [~, ~, V] = rsvd(mat_trans', num_of_pc);
    score = mat_trans' * V;
    disp_texp = ['preforming pca - done! ', num2str(s), '/', num2str(num_of_steps), ' times'];
    disp(disp_texp)
    [knn_idx, ~] = knnsearch(score, score, 'K', k_step + 1);
    disp_texp = ['calculating knn ', num2str(s), '/', num2str(num_of_steps), ' times'];
    disp(disp_texp)
    for cell = 1:size(mat_smooth, 2)
        mat_smooth(:, cell) = sum(raw_mat(:, knn_idx(cell, :)), 2);
    end
end

end
