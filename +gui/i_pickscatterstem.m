function [methodid] = i_pickscatterstem(defaultmethod, parentfig)
if nargin < 2, parentfig = []; end
if nargin < 1
    defaultmethod = 'Scatter+Stem';
end

methodid = [];

methodx = gui.myQuestdlg(parentfig, 'Plot type:', '', {'Scatter', 'Stem', ...
    'Scatter+Stem'}, defaultmethod);

if isempty(methodx), return; end

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