function [methodid] = i_pickscatterstem(defaultmethod, FigureHandle)
if nargin < 2, FigureHandle = []; end
if nargin < 1
    defaultmethod = 'Scatter+Stem';
end
methodx = gui.myQuestdlg(FigureHandle, 'Plot type:', '', {'Scatter', 'Stem', ...
    'Scatter+Stem'}, defaultmethod);
switch methodx
    case 'Scatter'
        methodid = 2;
    case 'Stem'
        methodid = 1;
    case 'Scatter+Stem'
        methodid = 5;
    otherwise
        methodid = [];
end
end