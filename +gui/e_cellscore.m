function [cs]=e_cellscore(sce,posg)

    answer = questdlg('Which method?',...
    'Select Method', 'AddModuleScore/Seurat', ...
    'UCell [PMID: 34285779]','AddModuleScore/Seurat');
    switch answer
        case 'AddModuleScore/Seurat'
            fw=gui.gui_waitbar;
            try
                [cs]=sc_cellscore(sce.X,sce.g,posg);
            catch ME
                gui.gui_waitbar(fw,true);
                errordlg(ME.message);                
            return;
            end
            gui.gui_waitbar(fw);
        case 'UCell [PMID: 34285779]'
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
    fprintf('\n=============\n');
    for k=1:length(posg)
        fprintf('%s\n',posg(k));
    end
    fprintf('=============\n');
end