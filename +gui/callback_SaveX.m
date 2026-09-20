function [OKPressed] = callback_SaveX(src, ~)

if isa(src, 'matlab.ui.Figure')
    FigureHandle = src;
    sce = guidata(FigureHandle);
else
    [FigureHandle, sce] = gui.gui_getfigsce(src);
end

% Capture the look now rather than trusting it to be current. Marker and
% colormap are recorded when they are picked, but the camera angle also
% moves by dragging the plot and through the many early-return branches of
% the 2D/3D switch, none of which announce themselves. Reading it here
% means what gets saved is what is on screen.
gui.i_storedisplay(src);

[OKPressed] = gui.sc_savescedlg(sce, FigureHandle);
end
