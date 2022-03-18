function [OKPressed]=callback_SaveX(src,~)

OKPressed=false;

FigureHandle=src.Parent.Parent;
sce=guidata(FigureHandle);

if ~(ismcc || isdeployed)
    answer = questdlg('Export & save data to:','',...
        'Workspace','MAT file','Seurat/RDS file','Workspace');
else
%    answer = questdlg('Export & save data to:','',...
%        'MAT file','Seurat/RDS file','MAT file');
    answer = questdlg('Export & save data to:','',...
        'Workspace','MAT file','Seurat/RDS file','Workspace');
end

switch answer
    case 'Workspace'
        labels = {'Save SCE to variable named:',...
            'Save SCE.X to variable named:',...
            'Save SCE.g to variable named:',...
            'Save SCE.S to variable named:'}; 
        vars = {'sce','X','g','s'};
        values = {sce,sce.X,sce.g,sce.s};
        [~,OKPressed]=export2wsdlg(labels,vars,values,...
            'Save Data to Workspace',...
            logical([1 0 0 0]),{@smhelp});
    case 'MAT file'
        [file, path] = uiputfile({'*.mat';'*.*'},'Save as');
        if isequal(file,0) || isequal(path,0)
           return;
        else
           filename=fullfile(path,file);
           fw=gui.gui_waitbar;
           save(filename,'sce','-v7.3');
           gui.gui_waitbar(fw);
           OKPressed=true;
        end
    case 'Seurat/RDS file'
        [file, path] = uiputfile({'*.rds';'*.*'},'Save as');
        if isequal(file,0) || isequal(path,0)
           return;
        else
           filename=fullfile(path,file);
           fw=gui.gui_waitbar;
           sc_sce2rds(sce,filename);
           gui.gui_waitbar(fw);
           disp("A<-readRDS(""input.rds"")");
           OKPressed=true;
        end
    otherwise
        return;
end         
    function smhelp
        helpdlg({'Select one or both check boxes.',...
                 'Change the variable names, if desired,',...
                 'and then click OK.'});
    end
end
