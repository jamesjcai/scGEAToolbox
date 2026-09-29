function callback_CellQCViolin(src, ~)
%CALLBACK_CELLQCVIOLIN Violin plots of the cell QC metrics, over all cells or by group.
%
%   Behind Plots > Cell QC Metrics in Violin Plots. Asks whether to pool all
%   cells or plot one violin per group, and for a grouping variable in the
%   second case.
%
%   See also GUI.I_QCVIOLIN, GUI.I_SELECT1CLASS.

[FigureHandle, sce] = gui.gui_getfigsce(src);

answer = gui.myQuestdlg(FigureHandle, ...
    'Plot the QC metrics over all cells, or one violin per group?', ...
    'Cell QC Metrics', {'All Cells', 'By Group', 'Cancel'}, 'By Group');
switch answer
    case 'All Cells'
        gui.i_qcviolin(sce.X, sce.g, FigureHandle);
    case 'By Group'
        [thisc, clabel] = gui.i_select1class(sce, false, [], [], FigureHandle);
        if isempty(thisc), return; end
        gui.i_qcviolin(sce.X, sce.g, FigureHandle, thisc, clabel);
    otherwise
        % Cancel or dismissed: nothing to plot.
end
end
