function [answer]=i_exporttable(T,needwait,TName,defname)
if nargin<4, defname=[]; end
if nargin<3, TName='T'; end
if nargin<2, needwait=false; end
    
if ~(ismcc || isdeployed)
    answer = questdlg('Export & save data to:','',...
        'Workspace','Text file','Excel file','Workspace');
else
    answer = questdlg('Export & save data to:','',...
        'Text file','Excel file','MAT file','Text file');
end
	
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
        if ~isempty(defname)
            [file, path] = uiputfile({'*.txt';'*.*'},'Save as',defname);
        else
            [file, path] = uiputfile({'*.txt';'*.*'},'Save as');
        end
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
        if ~isempty(defname)
            [file, path] = uiputfile({'*.xlsx';'*.xls';'*.*'}, ...
                'Save as',defname);
        else            
            [file, path] = uiputfile({'*.xlsx';'*.xls';'*.*'}, ...
                'Save as');
        end
        if isequal(file,0) || isequal(path,0), return; end
        
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
    case 'MAT file'
        [file, path] = uiputfile({'*.mat';'*.*'},'Save as');
        if isequal(file,0) || isequal(path,0), return; end
       filename=fullfile(path,file);
       %fw=gui.gui_waitbar;
       save(filename,'T');
       %gui.gui_waitbar(fw);
    otherwise
        return;
end
end
