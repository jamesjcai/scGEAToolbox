function [answer, id] = i_pickscoremethod(methodid, parentfig)

if nargin < 2, parentfig = []; end

if nargin < 1 || isempty(methodid)
    methodid = 0;
end

if methodid == 1
    answer = 'UCell [PMID:34285779]';
    id = 1;
elseif methodid == 2
    answer = 'AddModuleScore/Seurat';
    id = 2;
else
    answer = gui.myQuestdlg(parentfig, 'Select algorithm:', ...
        '', ...
        {'AddModuleScore/Seurat','UCell [PMID:34285779]'},...
        'AddModuleScore/Seurat');

    switch answer
        case 'UCell [PMID:34285779]'
            id = 1;
%           answer = 'UCell [PMID:34285779]';
        case 'AddModuleScore/Seurat'
            id = 2;
        otherwise
            id = [];
    end
end
