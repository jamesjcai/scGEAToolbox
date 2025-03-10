function [OKPressed] = callback_SaveX(src, ~)

% OKPressed = false;
if isa(src, 'matlab.ui.Figure')
    FigureHandle = src;
    sce = guidata(FigureHandle);
else
    [FigureHandle, sce, isui] = gui.gui_getfigsce(src);
end


[OKPressed] = gui.sc_savescedlg(sce);
end
