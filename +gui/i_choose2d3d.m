function [ndim] = i_choose2d3d(FigureHandle)

if nargin<1, FigureHandle = []; end
ndim = [];
answer3 = gui.myQuestdlg(FigureHandle, '3D or 2D?', '', ...
    {'3D', '2D'}, '3D');
switch answer3
    case '3D'
        ndim = 3;
    case '2D'
        ndim = 2;
    otherwise
        return;
end