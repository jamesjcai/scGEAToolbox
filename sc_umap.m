function s = sc_umap(X, ndim, donorm, dolog1p)
%UMAP embedding of cells
% s=sc_umap(X,3);

%see also: SC_TSNE, SC_PHATE
% s_phate=run.PHATE(X,3,true);
% s_umap=run.UMAP(X,3);

    if nargin < 4, dolog1p = true; end
    if nargin < 3, donorm = true; end
    if nargin < 2, ndim = 3; end
    % if ~issparse(X), X=sparse(X); end
    
    if donorm
        X = sc_norm(X, 'type', 'libsize');
        disp('Library-size normalization...done.')
    end
    if dolog1p
        X = log1p(X);
        disp('log1p transformation...done.')
    end

    done = false;
    try
        s = run.ml_UMAP(X, ndim);
        done = true;
    catch ME
        warning(ME.message);        
    end
    if ~done
        disp('Try older version Matlab UMAP (v4.4).');
        s = run.ml_UMAP(X, ndim, 44);
    end

end