function obj = embedcells(obj, methodtag, forced, usehvgs, ...
    ndim, numhvg, whitelist, showwaitbar)
if nargin < 8, showwaitbar = true; end
if nargin < 7, whitelist = []; end
if nargin < 6 || isempty(numhvg), numhvg = 2000; end
if nargin < 5 || isempty(ndim), ndim = 3; end
if nargin < 4 || isempty(usehvgs), usehvgs = true; end
if nargin < 3 || isempty(forced), forced = false; end
if nargin < 2, methodtag = 'tsne3d'; end

% USEHVGS takes a mode name as well as a logical, so a caller that only
% passes a choice through - GUI.I_GETHVGNUM's answer reaching here from the
% main app, say - does not have to assemble the gene list itself. See
% PKG.I_VALIDGENEMODE for the names.
[~, usehvgs, usemarkers, markerweight] = pkg.i_validgenemode(usehvgs);

if isempty(obj.s) || forced
    if isstring(methodtag) || ischar(methodtag)
        methodtag = lower(methodtag);
    end

    if usehvgs && size(obj.X, 1) > numhvg
        % The gene set is shared with PKG.E_CELLPCS; the ranking and its
        % rationale live in PKG.I_SELECTHVGS.
        keep = pkg.i_selecthvgs(obj.X, obj.g, numhvg);

        X = obj.X(keep, :);
        g = obj.g(keep);
    else
        X = obj.X;
        g = obj.g;
    end

    if usemarkers
        markerlist = pkg.i_getmarkerwhitelist(obj.g, obj.X);
        if isempty(whitelist)
            whitelist = markerlist;
        else
            whitelist = union(string(whitelist(:)), markerlist);
        end
    end

    nappended = 0;
    if ~isempty(whitelist)
        ngene = numel(g);
        [X, g] = pkg.i_appendgenes(X, g, whitelist, obj.X, obj.g);
        nappended = numel(g) - ngene;
        fprintf('EMBEDCELLS: %d additional whitelisted genes included.\n', ...
            nappended);
    end

    % Up-weighting the appended genes has to happen AFTER library-size
    % normalisation, so this normalises here and tells the embedder not to
    % do it again. Scaling raw counts instead does not work: it changes
    % every cell's total, which SC_NORM then divides back out unevenly.
    % Measured at w=5 on 100 genes, pre-normalisation scaling left the
    % weighted genes at 1.6x the baseline variance instead of the intended
    % 25x, and dragged the genes that were NOT weighted to 0.73x.
    %
    % PKG.I_APPENDGENES puts the added genes at the bottom, so they are the
    % last NAPPENDED rows.
    applyweight = markerweight > 1 && nappended > 0;
    if applyweight && ~any(strncmp(methodtag, {'tsne', 'umap'}, 4))
        % SC_PHATE and RUN.ML_METAVIZ normalise internally and expose no
        % flag to stop them, so there is nowhere to put the weights. Warn
        % rather than error: this runs per method inside a multi-method
        % loop, and aborting would lose the embeddings that do support it.
        warning('SingleCellExperiment:embedcells:weightNotApplied', ...
            ['%s normalises internally, so the marker weighting was not ' ...
            'applied to it. Its gene set is still the union.'], ...
            upper(methodtag));
        applyweight = false;
    end
    if applyweight
        X = log1p(sc_norm(X));
        X(isnan(X)) = 0;
        wrows = (numel(g) - nappended + 1):numel(g);
        X(wrows, :) = X(wrows, :) * markerweight;
        fprintf('EMBEDCELLS: those %d genes up-weighted %gx.\n', ...
            nappended, markerweight);
    end

    switch methodtag
        case {'tsne','tsne2d','tsne3d'}
            obj.s = sc_tsne(X, ndim, ~applyweight, ~applyweight);
        case {'umap','umap2d','umap3d'}
            obj.s = sc_umap(X, ndim, ~applyweight, ~applyweight);
        case {'phate','phate2d','phate3d'}
            obj.s = sc_phate(X, ndim);
        case {'metaviz','metaviz2d','metaviz3d'}
            obj.s = run.ml_metaviz(X, ndim, showwaitbar);
    end

    if contains(methodtag,'2d') || contains(methodtag,'3d')
        methoddimtag = methodtag;
    else
        % Bare method name ('tsne', 'umap', ...): append the dimension so the tag
        % matches the field names PKG.E_MAKEEMBEDSTRUCT defines. NDIM of 2 or 3
        % lands on an existing field ("tsne" with NDIM=3 becomes "tsne3d"); a
        % higher NDIM adds a new one, as before.
        methoddimtag = sprintf('%s%dd',methodtag, ndim);
    end

    obj.struct_cell_embeddings.(methoddimtag) = single(obj.s);
else
    disp('Use `sce = sce.embedcells(''tSNE'', true)` to overwrite existing SCE.S.');
end
end
