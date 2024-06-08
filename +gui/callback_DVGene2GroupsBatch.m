function callback_DVGene2GroupsBatch(src, ~)


FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

if ~gui.gui_showrefinfo('DV in Batch Mode'), return; end

if isscalar(unique(sce.c_cell_type_tx))
    warndlg('Only one cell type or cell type is undetermined.','');
    return;
end

warndlg('Function is under development.','');

return;
