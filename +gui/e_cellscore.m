function [cs]=e_cellscore(sce,posg,methodid)

cs=[];
if nargin<3 || isempty(methodid)
    answer = questdlg('Select algorithm:',...
    'Select Method', ...
    'UCell [PMID:34285779]','AddModuleScore/Seurat', ...
    'UCell [PMID:34285779]');
    switch answer
        case 'AddModuleScore/Seurat'
            methodid=2;
        case 'UCell [PMID:34285779]'
            methodid=1;            
        otherwise
            return;
    end
end

    fw=gui.gui_waitbar;
    try
        if methodid==1
            [cs]=sc_cellscore_ucell(sce.X,sce.g,posg);
        elseif methodid==2
            [cs]=sc_cellscore_admdl(sce.X,sce.g,posg);
        end                
    catch ME
        gui.gui_waitbar(fw,true);
        errordlg(ME.message);
        return;
    end
    gui.gui_waitbar(fw);



    posg=sort(posg);
    isexpressed=ismember(upper(posg),upper(sce.g));
    fprintf('\n=============\n%s\n-------------\n','Genes');
    for k=1:length(posg)
        if isexpressed(k)
           fprintf('%s\t*\n',posg(k));
        else
            fprintf('%s\t\n',posg(k));
        end
    end
    fprintf('=============\n*expressed genes (n=%d)\n', ...
        sum(isexpressed));

end