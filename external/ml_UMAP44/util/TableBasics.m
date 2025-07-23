classdef TableBasics<handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

properties(Constant)
    FILE_FILTERS={'*.xls', 'Excel (*.xls)'; ...
                '*.txt', 'Tab delimited text (*.txt)';...
                '*.csv', 'Comma separated values (*.csv)'};
    DEBUG=false;
end
    methods(Static)
        function resultSet=Select(T, column, value, ...
                columnsInResultSet, operation, resultSetIsArray)
            if nargin<6
                resultSetIsArray=[];
                if nargin<5
                    operation=[];
                    if nargin<4
                        columnsInResultSet=[];
                        if nargin<3
                            value=[];
                            if nargin<2
                                column=[];
                                if nargin<1
                                    T=[];
                                end
                            end
                        end
                    end
                end
            end
            if ~isa(T, 'table')
                error(['First argument to TableBasics.Select ' ...
                    'is a %s instead of the required ' ...
                    'matlab table object '], class(T));
            end
            [R, C]=size(T);
            if isempty(resultSetIsArray)
                resultSetIsArray='matrix';
            else
                resultSetIsArray=lower(resultSetIsArray);
            end
            allowedResultSetTypes={'matrix', 'table', 'rows'};
            if ~contains(allowedResultSetTypes, resultSetIsArray)
                error('resultSetType "%s" NOT allowed, allowed are %s', ...
                    resultSetIsArray, StringArray.toString( ...
                    allowedResultSetTypes, ',', true));
            end
            if isempty(value) && ~isempty(column) ...
                    && isempty(operation) && isempty(columnsInResultSet)
                if strcmp(resultSetIsArray, 'matrix')
                    resultSet=T{:, column};
                elseif strcmp(resultSetIsArray, 'table')
                    resultSet=T(:,column);
                else
                    resultSet=1:R;
                end
                return;
            end
            if isempty(column) || strcmp('*', column)
                if ~isempty(columnsInResultSet)
                    resultSet=T(:, columnsInResultSet);
                else
                    warning('Nothing to do on %dx%d table cells', R, C);
                    resultSet=T;
                end
                return;
            end
            if ~isempty(operation)
                operation=lower(operation);
                allowedOperations={'==', '>', '>=', '<', '<=', ...
                    'isnan', 'isempty', '~isnan', '~isempty'};
                if ~contains(allowedOperations, operation)
                    error('operation "%s" NOT allowed, allowed are %s', ...
                        operation, StringArray.toString( ...
                        allowedOperations, ',', true));
                end
            end
            isString=ischar(value);
            if isString
                left=lower(T{:, column});
                if strcmpi(operation, 'isempty')
                    rows=strcmp(left, '');
                elseif strcmpi(operation, '~isempty')
                    rows=~strcmp(left, '');
                else
                    doEndsWith=startsWith(value, '*');
                    doStartsWith=endsWith(value, '*');
                    if ~doStartsWith && ~doEndsWith
                        rows=strcmpi(left, value);
                    else
                        right=lower(strrep(value, '*',''));
                        if doStartsWith && doEndsWith
                            rows=contains(left, right);
                        elseif doStartsWith
                            rows=startsWith(left, right);
                        else
                            rows=endsWith(left, right);
                        end
                    end
                end
            else
                switch operation
                    case '>'
                        rows=T{:, column}>value;
                    case '<'
                        rows=T{:, column}>value;
                    case '>='
                        rows=T{:, column}>value;
                    case '<='
                        rows=T{:, column}>value;
                    case 'isnan'
                        rows=isnan(T{:, column}) | isempty(T{:, column});
                    case 'isempty'
                        rows=isnan(T{:, column}) | isempty(T{:, column});
                    case '~isnan'
                        rows=~isnan(T{:, column}) & ~isempty(T{:, column});
                    case '~isempty'
                        rows=~isnan(T{:, column}) & ~isempty(T{:, column});
                    otherwise
                        rows=T{:, column}==value;
                end
            end
            if isempty(columnsInResultSet)
                columnsInResultSet=column;
            elseif isequal(columnsInResultSet, '*')
                columnsInResultSet=1:C;
            end
            
            if strcmp(resultSetIsArray, 'table')
                resultSet=T(rows, columnsInResultSet);
            elseif strcmp(resultSetIsArray, 'matrix')
                try
                    resultSet=T{rows, columnsInResultSet};
                catch ex
                    warning(ex.identifier, 'Returning table because "%s"', ex.message);
                    resultSet=T(rows, columnsInResultSet);
                end
            else
                resultSet=rows;
            end
        end

        function filters=Filters(props, prop)
            dflt=props.get(prop, 'match_table.xls');
            filters=TableBasics.FILE_FILTERS;
            if endsWith(lower(dflt), '.csv')
                filters={filters{3, 1}, filters{3, 2};...
                    filters{1,1}, filters{1,2};...
                    filters{2,1}, filters{2,2}};
            elseif endsWith(lower(dflt), '.txt')
                filters={filters{2,1}, filters{2,2};...
                    filters{1,1}, filters{1,2};...
                    filters{3,1}, filters{3,2}};
            end
        end

        function file=UiPutFile(propFldr, propFile, props, dfltFile, rootFolder)
            if nargin<5
                rootFolder=File.Documents;
                if nargin<4
                    dfltFile='match_table.xls';
                    if nargin<3
                        props=BasicMap.Global;
                        if nargin<2
                            propFile=ClassificationTable.PROP_FILE;
                            if nargin<1
                                propFldr=ClassificationTable.PROP_FOLDER;
                            end
                        end
                    end
                end
            end
            dflt=props.get(propFile, dfltFile);
            filters={'*.xls', 'Excel (*.xls)'; ...
                '*.txt', 'Tab delimited text (*.txt)';...
                '*.csv', 'Comma separated values (*.csv)'};
            if endsWith(lower(dflt), '.csv')
                filters={filters{3, 1}, filters{3, 2};...
                    filters{1,1}, filters{1,2};...
                    filters{2,1}, filters{2,2}};
            elseif endsWith(lower(dflt), '.txt')
                filters={filters{2,1}, filters{2,2};...
                    filters{1,1}, filters{1,2};...
                    filters{3,1}, filters{3,2}};
            end
            [fldr, file]=uiPutFile(rootFolder, dflt, ...
                props, propFldr, 'Save as xls, txt or csv',...
                true, filters);
            if isempty(fldr)
                file=[];
                return;
            end
            props.set(propFile, file);
            file=fullfile(fldr, file);
        end

        function file=UiGetFile(propFldr, propFile, props, dfltFile, rootFolder, usePutFileIf1stTime)
            if nargin<6
                usePutFileIf1stTime=true;
                if nargin<5
                    rootFolder=File.Documents;
                    if nargin<4
                        dfltFile='5_classifiers.xls';
                        if nargin<3
                            props=BasicMap.Global;
                            if nargin<2
                                propFile=ClassificationTable.PROP_FILE;
                                if nargin<1
                                    propFldr=ClassificationTable.PROP_FOLDER;
                                end
                            end
                        end
                    end
                end
            end
            if usePutFileIf1stTime
                fl=props.get(propFldr);
                if isempty(fl)
                    file=TableBasics.UiPutFile(propFldr, propFile, props, dfltFile, rootFolder);
                    return;
                end
            end
            dflt=props.get(propFile, dfltFile);
            filters={'*.xls', 'Excel (*.xls)'; ...
                '*.txt', 'Tab delimited text (*.txt)';...
                '*.csv', 'Comma separated values (*.csv)'};
            if endsWith(lower(dflt), '.csv')
                filters={filters{3, 1}, filters{3, 2};...
                    filters{1,1}, filters{1,2};...
                    filters{2,1}, filters{2,2}};
            elseif endsWith(lower(dflt), '.txt')
                filters={filters{2,1}, filters{2,2};...
                    filters{1,1}, filters{1,2};...
                    filters{3,1}, filters{3,2}};
            end
            file=uiGetFile(filters, rootFolder, ...
                'Select xls, txt or csv', props, propFldr);
            if ~isempty(file)
                [~, f, e]=fileparts(file);
                props.set(propFile, [f e]);
            end
        end

        function html=Html(T, show)
            if nargin<2
                show=true;
            end
            [R,C]=size(T);
            if isempty(T.Properties.VariableDescriptions)
                H=T.Properties.VariableNames;
            else
                H=T.Properties.VariableDescriptions;
            end
            sb=java.lang.StringBuilder(R*C*15);
            sb.append("<table border='2' cellpadding='0' cellspacing='0'><thead><tr>");
            for c=1:C
                sb.append("<th>&nbsp;");
                sb.append(H{c});
                sb.append("&nbsp;</th>");
            end
            sb.append("</tr></thead>");
            if R>0
                conv=cell(1,C);
                for c=1:C
                    value=T{1,c};
                    if isnumeric(value)
                        conv{c}=@nCell;
                    else
                        conv{c}=@sCell;
                    end
                end
                for r=1:R
                    sb.append("<tr>");
                    for c=1:C
                        value=T{r,c};
                        feval(conv{c}, value);
                    end
                    sb.append("</tr>");
                end
            end
            sb.append("</table>");
            html=char(sb.toString);
            if show
                Html.BrowseString(Html.Wrap(html))
            end

            function nCell(x)
                if isnan(x)
                    sb.append('<td>&nbsp;N/A</td>');
                else
                    sb.append('<td align="right">');
                    sb.append(String.encodeRounded(x, 3, true));
                    sb.append("&nbsp;</td>");
                end
            end

            function sCell(x)
                x=x{1};
                if startsWith(x, "<html>")
                    sb.append("<td align='center'>");
                    sb.append(Html.remove(x));
                else
                    sb.append("<td>&nbsp;");
                    sb.append(Html.Remove(x));
                end
                sb.append("&nbsp;</td>");
            end
        end

        function data=NoIntsOrNans(T, cols)
            C=length(cols);
            data=[];
            for i=1:C
                idx=cols(i);
                numbs=double(T{:,idx});
                numbs(isnan(numbs))=0;
                data=[data numbs];
            end
        end

        function html=Summarize(T, names, fcns, nameColspan, cols, conv, ...
                trStarts, colSpan1)
            if nargin<8
                colSpan1=nameColspan;
                if nargin<7
                    trStarts={"<tr>", "<tr>"};
                end
            end
            sb=java.lang.StringBuilder;
            nFcns=length(fcns);
            C=length(cols);
            sb.append("<tr><td colspan='");
            sb.append(nameColspan+C);
            sb.append("' align='center' bgcolor'#AAAAAA'></td></tr>");
            data=TableBasics.NoIntsOrNans(T, cols);
            for f=1:nFcns
                nums=feval(fcns{f}, data);
                sb.append(trStarts{f});
                sb.append("<td align='right' colspan='");
                sb.append(num2str(colSpan1));
                sb.append("'><b>");
                sb.append(names{f});
                sb.append("</b>&nbsp;</td>");
                for c=1:C
                    feval(conv, nums(c), cols(c), sb);
                end
                sb.append('</tr>');
            end
            html=char(sb.toString);
        end

        function s=ToString(T, nameEqualsValueFmt)
            if nargin<2
                nameEqualsValueFmt=true;
            end
            if nameEqualsValueFmt
                [R,C]=size(T);
                s='';
                v=T.Properties.VariableNames;
                for r=1:R
                    for c=1:C
                        value=T{r,c};
                        if isnumeric(value)
                            s=[s v{c} '=' num2str(value)];
                        else
                            s=[s v{c} '-' value];
                        end
                        if c<C
                            s=[s ', '];
                        end
                    end
                    s=[s newline]; %#ok<*AGROW>
                end
            else
                s=formattedDisplayText(T);
            end
        end
    
        function [modelColumnIdxs, modelColumnNames, columnNames, ...
                columnNamesReordered, reordered]...
                =Match(modelColumnNames, columnNames, speakToUser)
            if nargin<3
                speakToUser=true;
            end
            columnNamesReordered=columnNames;
            [same, reordered]=StringArray.AreSameOrEmpty(...
                modelColumnNames, columnNames);
            if ~same
                [allFound, modelNameIdxs]=StringArray.Find(...
                    modelColumnNames, columnNames, true);
                if allFound
                    columnNames=columnNames(modelNameIdxs);
                    reordered=true;
                end
            end
            if ~same && ~reordered
                columnNames=StringArray.RemoveStartingAt(columnNames, ':');
                columnNamesReordered=columnNames;
                [same, reordered]=StringArray.AreSameOrEmpty(...
                    modelColumnNames, columnNames);
                if ~same
                    [allFound, modelNameIdxs]=StringArray.Find(...
                        modelColumnNames, columnNames, true);
                    if allFound
                        columnNames=columnNames(modelNameIdxs);
                        reordered=true;
                        speakWarning
                    else
                        mf=StringArray.RemoveStartingAt( ...
                            modelColumnNames, ':');
                        [same, reordered]=StringArray.AreSameOrEmpty(...
                            mf, columnNames);
                        if ~same
                            [allFound, modelNameIdxs]=StringArray.Find(...
                                mf, columnNames, true);
                            if allFound
                                columnNames=columnNames(modelNameIdxs);
                                modelColumnNames=mf;
                                reordered=true;
                                speakWarning
                            end
                        end
                    end
                end
            end
            if ~same
                if ~reordered
                    if speakToUser
                        msgError(Html.Wrap(['Test set <b>MUST'...
                            '</b> have <u>ALL</u> of the ' ...
                            '<b><font color="red">training ' ...
                            'set''s columns</font></b>:'  ...
                            Html.To2Lists(modelColumnNames, columnNames, ...
                            'ol', 'Training set', 'Test set', ...
                            2)]), 8, 'south west');
                    end
                    modelColumnNames=[];
                    modelColumnIdxs=[];
                else
                    [~, modelColumnIdxs]=StringArray.Find(...
                        modelColumnNames, columnNamesReordered, true);
                end
            end

            function speakWarning
                if speakToUser
                    msgWarning(Html.Wrap(['Incomplete ' ...
                        'matches of column names...<br><br>' ...
                        '<i>It is okay to continue ' ...
                        '</i> ... but MLP works<br>' ...
                        'best when data matches ' ...
                        'more exactly.']), 15, ...
                        'south east', 'Mismatch');
                end
            end
        end

        function [data, columnNames]=Reorder(modelColumnNames, data, columnNames)
            [~,b]=StringArray.Find(modelColumnNames, columnNames);
            tempData=data(:, b(end));
            tempColumn=columnNames{b(end)};
            data(:,end)=data(:,b(1));
            columnNames{end}=columnNames{b(1)};
            data(:,1)=tempData;
            columnNames{1}=tempColumn;

        end
    end
end