function [X] = sc_transform(X, varargin)

% https://www.biorxiv.org/content/10.1101/2021.06.24.449781v1.full
% acosh transformation based on the delta method
% shifted logarithm (log(x + c)) with a pseudo-count c, so that it approximates the acosh transformation
% randomized quantile and Pearson residuals

p = inputParser;
defaultType = 'PearsonResiduals';
validTypes = {'PearsonResiduals', 'kNNSmoothing', 'SCTransform', ...
    'FreemanTukey', 'csndm', 'SCT'};
checkType = @(x) any(validatestring(x, validTypes));

addRequired(p, 'X', @isnumeric);
addOptional(p, 'type', defaultType, checkType)
parse(p, X, varargin{:})


switch lower(p.Results.type)
    case 'pearsonresiduals'
        % analytic Pearson residuals
        % https://doi.org/10.1101/2020.12.01.405886
        % https://gist.github.com/hypercompetent/51a3c428745e1c06d826d76c3671797c

        u = (sum(X, 2) * sum(X, 1)) ./ sum(X(:));
        s = sqrt(u+(u.^2)./100);
        X = (X - u) ./ s;
        X(isnan(X)) = 0;
        n = size(X, 2);
        % clip to sqrt(n);
        sn = sqrt(n);
        X(X > sn) = sn;
        X(X < -sn) = -sn;
    case 'knnsmoothing'
        % K-nearest neighbor smoothing for high-throughput single-cell RNA-Seq data
        % https://doi.org/10.1101/217737
        %
        X = knn_smooth(X, 5, 10);

    case 'glmpca'

    case 'normalization_sqrt'
        % https://twitter.com/hippopedoid/status/1337028817219620864?s=20
    case 'sctransform'
        % sctransform: Variance Stabilizing Transformations for Single Cell UMI Data
        % Hafemeister & Satija 2019
        % https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
        [X] = run.r_SeuratSctransform(X);
    case 'csndm'
        [X] = run.mt_csndm_trans(X);
        %     case 'bigscale'
        %         pth=fullfile(pw1,'thirdparty/bigSCale');
        %         addpath(pth);
        %         % model=1. Log(x), then each row (gene) normalized between [-5:5]
        %         [X]=transform_bigscale(X);
    case 'sct'
        % sc_sct
        [X] = run.r_SeuratSctransform(X, string(1:size(X, 1)));
    case 'freemantukey'
        % https://github.com/flo-compbio/monet/blob/master/monet/util/expression.py
        % Applies the Freeman-Tukey transformation to stabilize variance."
        % https://www.biorxiv.org/content/10.1101/2020.06.08.140673v2.full
        % https://www.nature.com/articles/nmeth.2930
        X = sc_norm(X, 'type', 'deseq');
        X = sqrt(X) + sqrt(X+1);
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

%default
num_of_pc = 10;

% only want 1 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 1
    error('Requires at most 1 optional inputs');
end

% set defaults for optional inputs
optargs = num_of_pc;

% now put these defaults into the valuesToUse cell array,
% and overwrite the ones specified in varargin.
optargs = varargin;

% Place optional args in memorable variable names
[num_of_pc] = optargs{:};

mat_smooth = raw_mat;
num_of_steps = ceil(log2(k+1));
disp_text = ['number of steps: ', num2str(num_of_steps)];
disp(disp_text)
for s = 1:num_of_steps
    k_step = min(2^s-1, k);
    mat_tpm = median(sum(mat_smooth)) * bsxfun(@rdivide, mat_smooth, sum(mat_smooth));
    mat_trans = sqrt(mat_tpm) + sqrt(mat_tpm+1);
    [~, ~] = sort(sum(mat_smooth, 2), 'descend');
    disp_texp = ['preforming pca ', num2str(s), '/', num2str(num_of_steps), ' times'];
    disp(disp_texp)
    [~, ~, V] = rsvd(mat_trans', num_of_pc);
    score = mat_trans' * V;
    disp_texp = ['preforming pca - done! ', num2str(s), '/', num2str(num_of_steps), ' times'];
    disp(disp_texp)
    dist_matrix = squareform(pdist(score));
    disp_texp = ['calculating distance matrix ', num2str(s), '/', num2str(num_of_steps), ' times'];
    disp(disp_texp)
    for cell = 1:size(mat_smooth, 2)
        [~, c_sort_cell_idx] = sort(dist_matrix(cell, :));
        mat_smooth(:, cell) = sum(raw_mat(:, c_sort_cell_idx(1:k_step+1)), 2);
    end
end
disp done !
end
