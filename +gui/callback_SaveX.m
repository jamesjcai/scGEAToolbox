function [OKPressed] = callback_SaveX(src, ~)

    if isa(src, 'matlab.ui.Figure')
        FigureHandle = src;
        sce = guidata(FigureHandle);
    else
        [FigureHandle, sce] = gui.gui_getfigsce(src);
    end
    
    [OKPressed] = gui.sc_savescedlg(sce, FigureHandle);
end
