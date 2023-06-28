function [cs]=e_cellscore(sce,posg)

    answer = questdlg('Select algorithm:',...
    'Select Method', ...
    'UCell [PMID:34285779]','AddModuleScore/Seurat', ...
    'UCell [PMID:34285779]');
    switch answer
        case 'AddModuleScore/Seurat'
            fw=gui.gui_waitbar;
            try
                [cs]=sc_cellscore_admdl(sce.X,sce.g,posg);
            catch ME
                gui.gui_waitbar(fw,true);
                errordlg(ME.message);                
            return;
            end
            gui.gui_waitbar(fw);
        case 'UCell [PMID:34285779]'
            %[cs]=run.UCell(sce.X,sce.g,posg);
            fw=gui.gui_waitbar;
            try
                [cs]=sc_cellscore_ucell(sce.X,sce.g,posg);
            catch ME
                gui.gui_waitbar(fw,true);
                errordlg(ME.message);
            return;
            end
            gui.gui_waitbar(fw);
        otherwise
            return;
    end

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