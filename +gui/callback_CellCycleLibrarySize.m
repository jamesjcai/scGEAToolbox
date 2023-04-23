function callback_CellCycleLibrarySize(src,~)
    FigureHandle=src.Parent.Parent;    
    sce=guidata(FigureHandle);
    % if isempty(sce.c_cell_cycle_tx)
    %     warndlg('No cell cycle phase (SCE.C_CELL_CYCLE_TX is empty.)');
    %     return;
    % end
    % if length(unique(sce.c_cell_cycle_tx))==1
    %     if unique(sce.c_cell_cycle_tx)=="undetermined"
    %         warndlg('Cell cycle phase (SCE.C_CELL_CYCLE_TX) is undetermined.');
    %         return;
    %     end
    % end
    [sce]=gui.e_cellcycle(sce);    
    %pkg.i_violinplot(sum(sce.X),sce.c_cell_cycle_tx);
    [c,cL]=grp2idx(sce.c_cell_cycle_tx);
    gui.i_violinplot(sum(sce.X),cL(c),"",true,cL);    
    ylabel('Library Size (Sum of UMIs per cell)');
end
