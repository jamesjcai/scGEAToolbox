function sce = sc_fullembed(sce)
    % SC_FULLEMBED  Ensure all desired embeddings exist in single-cell object.
    %
    %   sce = sc_fullembed(sce) checks and renames legacy fields, then
    %   computes missing 2D and 3D embeddings for TSNE, UMAP, PHATE, and METAVIZ.

    % List of embedding methods and their dimensions
    methods = {'tsne', 'umap', 'phate', 'metaviz'};
    dims    = [2, 3];

    % Fix legacy fields (e.g., tsne → tsne3d)
    for i = 1:numel(methods)
        m = methods{i};
        fixLegacy(m, [m '3d']);
    end

    % Ensure all method/dimension combinations exist
    for i = 1:numel(methods)
        for d = dims
            newf = sprintf('%s%dd', methods{i}, d);
            embedIfMissing(newf, d);
        end
    end

    % ---- Subfunctions ----
    function fixLegacy(oldf, newf)
        s = sce.struct_cell_embeddings;
        if isfield(s, newf)
            return; % Already has the new field
        end
        if isfield(s, oldf) && ~isempty(s.(oldf)) && size(s.(oldf), 2) == 3
            s.(newf) = single(s.(oldf));  % rename 3D field
        end
        % Clean up empty or renamed old field
        if isfield(s, oldf) && isempty(s.(oldf))
            s = rmfield(s, oldf);
        elseif isfield(s, oldf) && strcmp(oldf, newf) == false
            s = rmfield(s, oldf);
        end
        sce.struct_cell_embeddings = s;
    end

    function embedIfMissing(fieldName, ndim)
        s = sce.struct_cell_embeddings;
        if ~isfield(s, fieldName) || isempty(s.(fieldName))
            sce = sce.embedcells(fieldName, true, true, ndim);
        end
    end
end

%{
function [sce]=sc_fullembed(sce)


in_fixfield('tsne','tsne3d');
in_fixfield('umap','umap3d');
in_fixfield('phate','phate3d');
in_fixfield('metaviz','metaviz3d');


in_embed('tsne2d',2);
in_embed('tsne3d',3);
in_embed('umap2d',2);
in_embed('umap3d',3);
in_embed('phate2d',2);
in_embed('phate3d',3);
% in_embed('metaviz2d',2);
% in_embed('metaviz3d',3);


    function in_embed(newf, ndim)
        if ~isfield(sce.struct_cell_embeddings,newf) || isempty(sce.struct_cell_embeddings.(newf))
            sce = sce.embedcells(newf, true, true, ndim);
        end   
    end

    function in_fixfield(oldf,newf)
        if ~isfield(sce.struct_cell_embeddings,newf) && isfield(sce.struct_cell_embeddings,oldf)
            if ~isempty(sce.struct_cell_embeddings.(oldf))
                if size(sce.struct_cell_embeddings.(oldf),2) == 3
                    sce.struct_cell_embeddings.(newf) = single(sce.struct_cell_embeddings.(oldf));
                    sce.struct_cell_embeddings = rmfield(sce.struct_cell_embeddings,oldf);
                end
            end
        end
        if isfield(sce.struct_cell_embeddings, oldf)
            if isempty(sce.struct_cell_embeddings.(oldf))
                sce.struct_cell_embeddings = rmfield(sce.struct_cell_embeddings,oldf);
            end
        end
    end

end

%}