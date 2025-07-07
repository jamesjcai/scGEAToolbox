function [ndim] = i_choose2d3dnmore(parentfig)

if nargin<1, parentfig = []; end
ndim = [];
answer = gui.myQuestdlg(parentfig, '3D or 2D?', '', ...
    {'3D', '2D', 'Other...'}, '3D');

if isempty(answer), return; end
switch answer
    case '3D'
        ndim = 3;
    case '2D'
        ndim = 2;
    case 'Other...'
        ndim = gui.i_inputnumk(20, 2, 100, ...
            'Number of Dimensions', parentfig);
    otherwise
        return;
end