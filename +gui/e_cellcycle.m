function [sce] = e_cellcycle(sce, parentfig)

if nargin<2, parentfig = []; end

needestimate = false;
if isempty(sce.c_cell_cycle_tx) || ...
        all(strcmp(unique(sce.c_cell_cycle_tx), "undetermined"))
    needestimate = true;
else
    answer1 = gui.myQuestdlg(parentfig, ...
        ['Use existing cell cycle estimation or ' ...
        're-compute new estimation?'], ...
        '', {'Use existing', 'Re-compute', 'Cancel'}, 'Use existing');
    switch answer1
        case 'Re-compute'
            needestimate = true;
        case 'Cancel'
            return;
    end
end
if needestimate
    fw = gui.gui_waitbar;
    sce = sce.estimatecellcycle(true, 1);
    gui.gui_waitbar(fw);
end
end
