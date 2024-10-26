function s = sc_tsne(X, ndim, donorm, dolog1p)
%tSNE embedding of cells

%see also: RUN.MT_PHATE, RUN.MT_UMAP
% s_phate=run.PHATE(X,3,true);
% s_umap=run.UMAP(X,3);
narginchk(1, 4)
% validateattributes(ndim, {'numeric'}, ...
%     {'scalar', 'integer', '>=', 2, '<=', 3});

if nargin < 2, ndim = 3; end
if nargin < 3, donorm = true; end
if nargin < 4, dolog1p = true; end

% if nargin<5, bygene=false; end   % when BYGENE=true, the matrix X will be transposed and the output will be tSNE for genes rather than cells.
% if nargin<6, genelist=[]; end

pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '+run', 'external', 'ml_PHATE');
if ~(ismcc || isdeployed), addpath(pth); end

%if bygene, X=X.'; end
if donorm
    X = sc_norm(X);
    disp('Library-size normalization...done.')
end

X(isnan(X)) = 0; % empty cells will be retained to keep spots

if dolog1p
    X = log1p(X);
    disp('log1p transformation...done.')
end

if issparse(X), X = full(X); end
[ngenes] = size(X, 1);

% The following transpose is necessary to make the input dim right.
data = X.';
%if ncells>500
%	data = svdpca(data, 50, 'random');
%end
    if ngenes > 500
        if ndim < 10
            s = tsne(data, 'NumDimensions', ndim, ...
                'Algorithm', 'barneshut', 'NumPCAComponents', 50, ...
                'Standardize', false);
        else
            s = tsne(data, 'NumDimensions', ndim, ...
                'Algorithm', 'barneshut', 'NumPCAComponents', 5*ndim, ...
                'Standardize', false);
        end
    else
        s = tsne(data, 'NumDimensions', ndim, ...
            'Algorithm', 'exact', 'NumPCAComponents', 0, ...
            'Standardize', false);
    end
end
