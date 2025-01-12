function [OKPressed] = callback_SaveX(src, ~)

% OKPressed = false;
if isa(src, 'matlab.ui.Figure')
    FigureHandle = src;
else
    FigureHandle = src.Parent.Parent;
end
sce = guidata(FigureHandle);

[OKPressed] = gui.sc_savescedlg(sce);
end
