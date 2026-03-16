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
elseif methodid == 3
    answer = 'AUCell (AUC recovery)';
    id = 3;
else
    answer = gui.myQuestdlg(parentfig, 'Select algorithm:', ...
        '', ...
        {'AddModuleScore/Seurat','UCell [PMID:34285779]','AUCell (AUC recovery)'},...
        'AddModuleScore/Seurat');

    switch answer
        case 'UCell [PMID:34285779]'
            id = 1;
        case 'AddModuleScore/Seurat'
            id = 2;
        case 'AUCell (AUC recovery)'
            id = 3;
        otherwise
            id = [];
    end
end
