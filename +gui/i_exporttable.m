function [answer, filename] = i_exporttable(T, needwait, TName, ...
    deffilename, outtype, sheetname, parentfig)

    if nargin < 7, parentfig = []; end
    if nargin < 6 || isempty(sheetname), sheetname = []; end
    if nargin < 5 || isempty(outtype), outtype = []; end
    if nargin < 4 || isempty(deffilename), deffilename = []; end
    if nargin < 3 || isempty(TName), TName = 'T'; end
    if nargin < 2 || isempty(needwait), needwait = true; end
    filename = [];
    
    if ~isempty(outtype)
        answer = outtype;
    else
        if ~(ismcc || isdeployed)
            answer = gui.myQuestdlg(parentfig, 'Export & save data to:', '', ...
                {'Workspace', 'Text file', 'Excel file'}, 'Workspace');
        else
            answer = gui.myQuestdlg(parentfig, 'Export & save data to:', '', ...
                {'Text file', 'Excel file', 'MAT file'}, 'Text file');
        end
    end
    
    if isstring(TName)
        TName = char(TName);
    end
    
    switch answer
        case 'Workspace'
            labels = {'Save to variable named:'};
            vars = {TName};
            values = {T};
    
            if gui.i_isuifig(parentfig)
                gui.myExport2wsdlg(labels, vars, values, ...
                    'Save Data to Workspace', [], parentfig);
            else
                if needwait
                    %disp('needwait')
                    waitfor(export2wsdlg(labels, vars, values));
                else
                    %disp('~needwait')
                    export2wsdlg(labels, vars, values);
                end
            end
    
    
        case 'Text file'
            if ~isempty(deffilename)
                [file, path] = uiputfile({'*.txt'; '*.*'}, 'Save as', deffilename);
            else
                [file, path] = uiputfile({'*.txt'; '*.*'}, 'Save as');
            end
            if isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure'), figure(parentfig); end
            if isequal(file, 0) || isequal(path, 0)
                return;
            else
                filename = fullfile(path, file);
                try
                    writetable(T, filename, 'Delimiter', '\t', 'WriteRowNames', true);
                catch
                    writematrix(T, filename, 'Delimiter', '\t');
                end
                pause(1);
                % if needwait
                %     gui.myHelpdlg(parentfig, sprintf('Result has been saved in %s', filename), '');
                % else
                    gui.myHelpdlg(parentfig, ...
                        sprintf('Result has been saved in %s', filename));
    %            end
            end
        case 'Excel file'
            if ~isempty(deffilename)
                [file, path] = uiputfile({'*.xlsx'; '*.xls'; '*.*'}, ...
                    'Save as', deffilename);
            else
                [file, path] = uiputfile({'*.xlsx'; '*.xls'; '*.*'}, ...
                    'Save as');
            end
            % if isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure'), figure(parentfig); end
            if isequal(file, 0) || isequal(path, 0), return; end
    
            filename = fullfile(path, file);
    
            %{
            variables = T.Properties.VariableNames;
            for k = 1:length(variables)
                xx = T.(variables{k});
                if isnumeric(xx) && any(isinf(xx))
                    xx(isinf(xx) & xx > 0) = 1e99;
                    xx(isinf(xx) & xx < 0) = -1e99;
                    T.(variables{k}) = xx;
                end
            end
            %}
            for varIdx = 1:width(T)
                if isnumeric(T{:, varIdx})
                    T{T{:, varIdx} == Inf, varIdx} = 1e99;
                    T{T{:, varIdx} == -Inf, varIdx} = -1e99;
                end
            end
            
            if isempty(sheetname)
                writetable(T, filename, 'FileType', 'spreadsheet', ...
                    'WriteRowNames', true);
            else
                writetable(T, filename, 'FileType', 'spreadsheet', ...
                    'WriteRowNames', true, 'Sheet', sheetname);
            end
            pause(1)
            if needwait
                gui.myQuestdlg(parentfig, ...
                    sprintf('Result has been saved in %s', filename),'', ...
                    {'OK'},'OK');
            end
        case 'MAT file'
            if ~isempty(deffilename)
                [file, path] = uiputfile({'*.mat'; '*.*'}, ...
                    'Save as', deffilename);
            else        
                [file, path] = uiputfile({'*.mat'; '*.*'}, 'Save as');
            end
            if isequal(file, 0) || isequal(path, 0), return; end
            filename = fullfile(path, file);
            save(filename, 'T');
        otherwise
            return;
    end
    
end