function callback_CellCycleLibrarySize(src,~)
    FigureHandle=src.Parent.Parent;    
    sce=guidata(FigureHandle);
    if isempty(sce.c_cell_cycle_tx)
        warndlg('No cell cycle phase (SCE.C_CELL_CYCLE_TX is empty)');
        return;
    end
    figure; 
    pkg.i_violinplot(sum(sce.X),sce.c_cell_cycle_tx);
    ylabel('Library Size (Sum of UMIs per cell)');
end
