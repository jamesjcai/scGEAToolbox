function callback_SaveX(src,~)
answer = questdlg('Export & save data to:','',...
    'Workspace','MAT file','Seurat/RDS file','Workspace');
        FigureHandle=src.Parent.Parent;
        sce=guidata(FigureHandle);
switch answer
    case 'Workspace'
        labels = {'Save SCE to variable named:',...
            'Save SCE.X to variable named:',...
            'Save SCE.g to variable named:',...
            'Save SCE.S to variable named:'}; 
        vars = {'sce','X','genelist','s'};
        values = {sce,sce.X,sce.g,sce.s};
        export2wsdlg(labels,vars,values,...
            'Save Data to Workspace',...
            logical([1 0 0 0]),{@smhelp});
    case 'MAT file'
        [filen, pathn] = uiputfile( ...
           {'*.mat';'*.*'},'Save as');
        filename=[pathn,filen];
        if ~(filename), return; end
        save(filename,'sce');
    case 'Seurat/RDS file'        
        [filen, pathn] = uiputfile( ...
           {'*.rds';'*.*'},'Save as');
        filename=[pathn,filen];
        if ~(filename), return; end
        fw=gui.gui_waitbar;
        [status]=sc_sce2rds(sce,filename);
        gui.gui_waitbar(fw);
    otherwise
        return;
end
         
    function smhelp
        helpdlg({'Select one or both check boxes.',...
                 'Change the variable names, if desired,',...
                 'and then click OK.'});
    end
end
