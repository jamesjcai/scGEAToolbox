function obj = embedcells(obj, methodtag, forced, usehvgs, ...
    ndim, numhvg, whitelist)
if nargin < 7, whitelist = []; end
if nargin < 6, numhvg = 2000; end
if nargin < 5, ndim = 3; end
if nargin < 4 || isempty(usehvgs), usehvgs = true; end
if nargin < 3 || isempty(forced), forced = false; end
if nargin < 2, methodtag = 'tsne3d'; end

% validTypes = {'tsne','umap','phate'};
% checkType = @(x) any(validatestring(x,validTypes));
if isempty(obj.s) || forced
    if isstring(methodtag) || ischar(methodtag)
        methodtag = lower(methodtag);
    end

    if usehvgs && size(obj.X, 1) > numhvg
        % disp('Identifying HVGs')
        [~, X, g] = sc_hvg(obj.X, obj.g, true, false, true, false, true);


        %size(X)
        %size(obj.X)


        idx = false(size(obj.g));
        idx(1:numhvg) = true;

        try
            pw1 = fileparts(mfilename('fullpath'));
            pth = fullfile(pw1, '..', '+run', 'thirdparty', 'alona_panglaodb','marker_hs.mat');
            load(pth,'Tw');
            markerg = unique(string(Tw.Var1));
            idx(ismember(upper(obj.g), upper(markerg))) = true;
            pth = fullfile(pw1, '..', '+run', 'thirdparty', 'alona_panglaodb','marker_mm.mat');
            load(pth,'Tw');
            markerg = unique(string(Tw.Var1));
            idx(ismember(upper(obj.g), upper(markerg))) = true;
            fprintf('EMBEDCELLS: %d additional marker genes included.\n', sum(idx(numhvg+1:end)));
        catch
            warning('EMBEDCELLS: marker genes not included.');
        end

        X = X(idx, :);
        g = g(idx);


    else
        X = obj.X;
        g = obj.g;
    end

    if ~isempty(whitelist)
        assert(all(ismember(whitelist, obj.g)));
        [~, idx] = setdiff(whitelist, g);
        if ~isempty(idx)
            [~, idxx] = ismember(whitelist, obj.g);
            Xresv = obj.X(idxx, :);
            X = [X; Xresv(idx, :)];
            g = [g; whitelist(idx)];
        end
    end

    switch methodtag
        case {'tsne','tsne2d','tsne3d'}
            obj.s = sc_tsne(X, ndim, true);
        case {'umap','umap2d','umap3d'}
            obj.s = sc_umap(X, ndim);
        case {'phate','phate2d','phate3d'}
            obj.s = sc_phate(X, ndim);
        case {'metaviz','metaviz2d','metaviz3d'}
            obj.s = run.mt_metaviz(X, ndim);
    end

    % if contains(methodtag,'2d') || contains(methodtag,'3d')
    %     methoddimtag = methodtag;
    % elseif ndim > 3
        methoddimtag = sprintf('%s%dd',methodtag, ndim);
    % end

    obj.struct_cell_embeddings.(methoddimtag) = single(obj.s);
    % disp('SCE.S added');
else    
    disp('Use `sce = sce.embedcells(''tSNE'', true)` to overwirte existing SCE.S.');
end
end
