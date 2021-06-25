function callback_CalculateGeneStats(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    [T]=sc_genestats(sce.X,sce.g)
    
answer = questdlg('Export & save data to:','',...
    'Workspace','Text file','Excel file','Workspace');
	
switch answer
    case 'Workspace'
        labels = {'Save T to variable named:'}; 
        vars = {'T'};
        values = {T};
        export2wsdlg(labels,vars,values);
    case 'Text file'
        [file, path] = uiputfile({'*.txt';'*.*'},'Save as');
        if isequal(file,0) || isequal(path,0)
           return;
        else			
           filename=fullfile(path,file);
		   writetable(T,filename,'Delimiter','\t');
           pause(1)
           helpdlg(sprintf('Result has been saved in %s',filename))
        end
    case 'Excel file'
        [file, path] = uiputfile({'*.xlsx';'*.xls';'*.*'},'Save as');
        if isequal(file,0) || isequal(path,0)
           return;
        else
           filename=fullfile(path,file);
		   writetable(T,filename,'FileType','spreadsheet');
           pause(1)
           helpdlg(sprintf('Result has been saved in %s',filename))
        end
    otherwise
        return;
end

end