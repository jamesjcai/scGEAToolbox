function i_exporttable(T,needwait,TName)
if nargin<3, TName='T'; end
if nargin<2, needwait=false; end
    
answer = questdlg('Export & save data to:','',...
    'Workspace','Text file','Excel file','Workspace');
	
switch answer
    case 'Workspace'
        labels = {'Save to variable named:'}; 
        vars = {TName};
        values = {T};
        if needwait
            waitfor(export2wsdlg(labels,vars,values));
        else
            export2wsdlg(labels,vars,values);
        end
    case 'Text file'
        [file, path] = uiputfile({'*.txt';'*.*'},'Save as');
        if isequal(file,0) || isequal(path,0)
           return;
        else			
           filename=fullfile(path,file);
		   writetable(T,filename,'Delimiter','\t');
           pause(1)
           if needwait
                waitfor(helpdlg(sprintf('Result has been saved in %s',filename),''));
           else
               helpdlg(sprintf('Result has been saved in %s',filename),'')
           end
        end
    case 'Excel file'
        [file, path] = uiputfile({'*.xlsx';'*.xls';'*.*'},'Save as');
        if isequal(file,0) || isequal(path,0)
           return;
        else
           filename=fullfile(path,file);
           
           variables=T.Properties.VariableNames;
           for k=1:length(variables)
               xx=T.(variables{k});
               if isnumeric(xx) && any(isinf(xx))
                   xx(isinf(xx)&xx>0)=1e99;
                   xx(isinf(xx)&xx<0)=-1e99;
                   T.(variables{k})=xx;
               end
           end
		   writetable(T,filename,'FileType','spreadsheet');
           pause(1)
           if needwait
                waitfor(helpdlg(sprintf('Result has been saved in %s',filename),''));
           else
               helpdlg(sprintf('Result has been saved in %s',filename),'')
           end
        end
    otherwise
        return;
end
end
