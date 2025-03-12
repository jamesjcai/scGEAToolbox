function [ndim] = i_choose2d3d(parentfig)

if nargin<1, parentfig = []; end
ndim = [];
answer3 = gui.myQuestdlg(parentfig, '3D or 2D?', '', ...
    {'3D', '2D'}, '3D');
switch answer3
    case '3D'
        ndim = 3;
    case '2D'
        ndim = 2;
    otherwise
        return;
end