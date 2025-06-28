function [OKPressed] = callback_SaveX(src, ~)

% OKPressed = false;
if isa(src, 'matlab.ui.Figure')
    FigureHandle = src;
    sce = guidata(FigureHandle);
elseif isa(src, 'matlab.apps.AppBase')    
    [FigureHandle, sce] = xui.gui_getfigsce(src);
else
    [FigureHandle, sce] = gui.gui_getfigsce(src);
end


[OKPressed] = gui.sc_savescedlg(sce, FigureHandle);
end
