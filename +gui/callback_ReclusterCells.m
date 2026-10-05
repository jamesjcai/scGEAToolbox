function [requirerefresh] = callback_ReclusterCells(src, ~, methodtag, sx)

requirerefresh = false;
if nargin < 4, sx = []; end
methodtag = lower(methodtag);

[FigureHandle, sce] = gui.gui_getfigsce(src);

usingold = false;
% ISFIELD, not a direct read: an SCE restored from a .mat file carries
% whatever STRUCT_CELL_CLUSTERINGS held when it was saved, so a method
% added since then is simply absent rather than empty. The length check
% skips a result left over from before cells were added.
hasold = isfield(sce.struct_cell_clusterings, methodtag) && ...
    numel(sce.struct_cell_clusterings.(methodtag)) == sce.NumCells;

% SEURAT is Seurat's own FindClusters result, stored by the R run behind
% Cluster > "Embed Cells with Seurat" and not applied there. It cannot be
% recomputed here, only applied.
if strcmp(methodtag, 'seurat')
    if ~hasold
        gui.myErrordlg(FigureHandle, ['There is no Seurat clustering to ' ...
            'use. Run Cluster > Embed Cells with Seurat (R) first; it stores ' ...
            'the clusters Seurat finds.'], '');
        return;
    end
    sce.c_cluster_id = sce.struct_cell_clusterings.seurat;
    gui.myGuidata(FigureHandle, sce, src);
    requirerefresh = true;
    return;
end

if hasold
    answer1 = gui.myQuestdlg(FigureHandle, sprintf(['Cells already have ' ...
        'a %s clustering. Use it, or compute a new one?'], ...
        i_methodlabel(methodtag)), '', ...
        {'Yes, use existing', 'No, re-compute', ...
        'Cancel'}, 'Yes, use existing');
    switch answer1
        case 'Yes, use existing'
            sce.c_cluster_id = sce.struct_cell_clusterings.(methodtag);
            usingold = true;
        case 'No, re-compute'
            usingold = false;
        otherwise
            % Cancel, or the dialog closed.
            return;
    end
end

if ~usingold
    % LOUVAINPC is the Seurat-like route: the resolution, not a preset k,
    % decides the number of clusters, so that is what it asks for first.
    % Asking for k instead tunes the resolution until k clusters come out.
    % Choosing by annotation takes the resolution SC_ANNOTATIONSTABILITY
    % recommends; its partition is the one Resolution= would give, so it is
    % used as it is rather than recomputed. A list, not GUI.MYQUESTDLG:
    % three choices plus the Cancel it appends is one more than it shows.
    resolution = [];
    if strcmp(methodtag, 'louvainpc')
        choices = ["Resolution (larger gives more clusters)", ...
            "Number of clusters", ...
            "Choose by annotation (resolution whose cell types are stable)"];
        [indx, tf] = gui.myListdlg(FigureHandle, choices, ...
            'Louvain on Principal Components', 1, false, false, [380, 150], ...
            'How should the number of clusters be set?');
        if tf ~= 1, return; end
        switch indx
            case 1
                resolution = gui.i_askresolution(FigureHandle);
                if isempty(resolution), return; end
                k = [];
            case 2
                k = i_asknumclusters(sce, FigureHandle);
                if isempty(k), return; end
            case 3
                S = gui.i_annotationsweep(src, FigureHandle, sce);
                if isempty(S), return; end
                best = S.T(S.T.Recommended, :);
                msg = sprintf(['Recommended resolution %.3g: %d clusters.' ...
                    '\n\nApply reclusters the cells at this resolution; cell ' ...
                    'types are not changed, so annotate the clusters ' ...
                    'afterwards. Not Now leaves the clusters as they are; ' ...
                    'the sweep is saved with the data, so running this ' ...
                    'again offers it without recomputing.'], ...
                    S.Resolution, best.NumClusters);
                if ~strcmp(gui.myQuestdlg(FigureHandle, msg, 'Apply Resolution', ...
                        {'Apply', 'Not Now'}, 'Apply'), 'Apply')
                    return;
                end
                sce.c_cluster_id = S.ClusterId;
                sce.struct_cell_clusterings.louvainpc = S.ClusterId;
                sce.setCellAttribute('annotation_stability', S.Stability);
                gui.myGuidata(FigureHandle, sce, src);
                requirerefresh = true;
                return;
            otherwise
                % myListdlg with allowmulti=false returns one of the three indices.
                return;
        end
    else
        k = i_asknumclusters(sce, FigureHandle);
        if isempty(k), return; end
    end

    fw = gui.myWaitbar(FigureHandle);
    try
        sce = sce.clustercells(k, methodtag, true, sx, Resolution=resolution);
    catch ME
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return
    end
    gui.myWaitbar(FigureHandle, fw);
end

gui.myGuidata(FigureHandle, sce, src);
requirerefresh = true;
end

function label = i_methodlabel(methodtag)
% The names the method lists show, not the STRUCT_CELL_CLUSTERINGS field.
switch methodtag
    case 'louvainpc'
        label = 'Louvain (principal components)';
    case 'louvain'
        label = 'Louvain (embedding S)';
    case 'kmeans'
        label = 'k-means';
    case 'snndpc'
        label = 'SnnDpc';
    case 'sc3'
        label = 'SC3';
    otherwise
        % Methods reached only programmatically keep their field name.
        label = upper(methodtag);
end
end

function k = i_asknumclusters(sce, FigureHandle)
defv = round(sce.NumCells/100, -1);
defva = min([2, round(sce.NumCells/100, -2), ...
    round(sce.NumCells/20, -1)]);
if defva == 0, defva = min([2, defv]); end
defvb = max([round(sce.NumCells/20, -2), ...
    round(sce.NumCells/20, -1)]);

if any([defv, defva, defvb]==0)
    defv=5; defva=1; defvb=100;
end
k = gui.i_inputnumk(defv, defva, defvb, ...
    'Enter the number of clusters', FigureHandle);
end
