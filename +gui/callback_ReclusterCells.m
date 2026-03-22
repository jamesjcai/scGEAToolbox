function [requirerefresh] = callback_ReclusterCells(src, ~, methodtag, sx)

requirerefresh = false;
if nargin < 4, sx = []; end
methodtag = lower(methodtag);

[FigureHandle, sce] = gui.gui_getfigsce(src);

usingold = false;
if ~isempty(sce.struct_cell_clusterings.(methodtag))
    answer1 = gui.myQuestdlg(FigureHandle, sprintf('Using existing %s clustering?', ...
        upper(methodtag)), '', ...
        {'Yes, use existing', 'No, re-compute', ...
        'Cancel'}, 'Yes, use existing');
    switch answer1
        case 'Yes, use existing'
            sce.c_cluster_id = sce.struct_cell_clusterings.(methodtag);
            usingold = true;
        case 'No, re-compute'
            usingold = false;
        case 'Cancel'
            return;
    end
end

if ~usingold
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

    if isempty(k), return; end
    fw = gui.myWaitbar(FigureHandle);
    try
        sce = sce.clustercells(k, methodtag, true, sx);
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
