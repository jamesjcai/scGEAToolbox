classdef ClassificationTable<handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

properties(Constant)
    PROP_FILE='ClassificationTable.File';
    PROP_FOLDER='ClassificationTable.Folder';
    COLUMNS={'Classifier', 'DataSet', 'TestCase', 'CellType', ...
        'Size', 'F1Score', 'InlierConcordance', 'Similarity', 'Concordance'};
    TYPES={'string', 'string', 'string', 'string', 'int64', 'double', 'double', 'double', 'double' };
    COLOR=[.2 .2 .5];
    RGB='\\color[rgb]{.2 .2 .5}';
    VERSION_COLORS='1.2';
    VERSION_CLASS='1.0';
end
    methods(Static)

        function [localFile, localFileSpell, uri, uriSpell, fileName, fileNameSpell]=ColorFiles()
            fileName=['colorsByName_v' num2str( ...
                ClassificationTable.VERSION_COLORS) '.properties'];
            fileNameSpell=['colorsByName_v' num2str( ...
                ClassificationTable.VERSION_COLORS) '.spell.properties'];
            localFile=File.Downloads('suh_pipelines', 'storage.googleapis.com',...
                'cytogenie.org', fileName);
            uri=['https://storage.googleapis.com/cytogenie.org/' fileName];
            localFileSpell=File.Downloads('suh_pipelines', 'storage.googleapis.com',...
                'cytogenie.org', fileNameSpell);
            uriSpell=['https://storage.googleapis.com/cytogenie.org/' fileNameSpell];
        end

        function [file, fileSpell]=DownloadColors(updateGlobals)
            if nargin<1
                updateGlobals=false;
            end
            [f, fSpell, uri, uriSpell]=ClassificationTable.ColorFiles;
            if exist(f, 'file')
                delete(f);
            end
            if exist(fSpell, 'file')
                delete(fSpell);
            end
            file=WebDownload.LocateUri(uri);
            if ~exist(file, 'file')
                file=[];
            end
            fileSpell=WebDownload.LocateUri(uriSpell);
            if ~exist(fileSpell, 'file')
                fileSpell=[];
            end
            if updateGlobals
                app=BasicMap.Global;
                f1=File.Documents('run_umap', 'examples', ...
                    'colorsByName.properties');
                f2=File.Documents('run_umap', 'examples', ...
                    'colorsByName.spell.properties');
                to=File.Downloads('suh_pipelines', 'bak', File.Time);
                File.mkDir(to);
                copyfile(f1, to);
                copyfile(f2, to);
                File.OpenFolderWindow(to, 'openFolder', true, false, ...
                    ' for color backups')
                copyfile(file, f1);
                copyfile(fileSpell, f2);
                app.setColorsByName(ColorsByName);
            end
        end

        function UploadColors()
            if ~isequal('/Users/swmeehan', File.Home) || ~ismac
                msgError('Must be run on Stephen Meehan''s mac')
                return;
            end
            app=BasicMap.Global;
            cbn=app.colorsByName;
            [~,~,~,~,fileName, fileNameSpell]=ClassificationTable.ColorFiles;
            File.UploadFileToGoogleCloud(cbn.file, fileName, 'cytogenie.org');
            File.UploadFileToGoogleCloud(cbn.spellFile, fileNameSpell, 'cytogenie.org');
        end

        function [localFile, uri, fileName]=DownloadClassFile(nClassifiers)
            fileName=[num2str(nClassifiers) '_classifiers.xls'];
            localFile=File.Downloads('suh_pipelines', ...
                'storage.googleapis.com',...
                'cytogenie.org', fileName);
            uri=['https://storage.googleapis.com/cytogenie.org/' fileName];
        end

        function [localFile, uri, fileName]=UploadClassFile(nClassifiers)
            fileName=[num2str(nClassifiers) '_classifiers.xls'];
            localFile=File.GoogleDrive('FlowJoBridge', 'Papers', ...
                'DataSets', fileName);
            if ~exist(localFile, 'file')
                msgError(Html.WrapHr(['<b>File not found!</b>' ...
                    Html.FileTree(localFile)]));
                assert(exist(localFile, 'file'));
            end
            uri=['https://storage.googleapis.com/cytogenie.org/' fileName];
        end

        function [file, preExisted]=DownloadClasses(nClassifiers, force)
            if nargin<2
                force=false;
                if nargin<1
                    nClassifiers=5;
                end
            end
            preExisted=false;
            [f, uri]=ClassificationTable.DownloadClassFile(nClassifiers);
            if exist(f, 'file')
                preExisted=true;
                if force
                    delete(f);
                else
                    file=f;
                    return;
                end
            end
            file=WebDownload.LocateUri(uri);
            if ~exist(file, 'file')
                file=[];
            end
        end

        function UploadClasses()
            go(5);
            go(6);
            function go(nClassifiers)
                if ~isequal('/Users/swmeehan', File.Home) || ~ismac
                    msgError('Must be run on Stephen Meehan''s mac')
                    return;
                end
                [from, ~, fileName]=ClassificationTable.UploadClassFile(nClassifiers);
                tn=tempname;
                to=fullfile(fileparts(tn), fileName);
                copyfile(from, to);
                fprintf('Copy of uploaded %s is saved in temporary folder %s\n', ...
                    fileName, fileparts(tn))
                File.UploadFileToGoogleCloud(to, fileName, 'cytogenie.org');
            end
        end

        function T=New(rowSize)
            size=[rowSize, length(ClassificationTable.COLUMNS)];
            T=table('Size', size, ...
                'VariableType', ClassificationTable.TYPES, ...
                'VariableNames', ClassificationTable.COLUMNS);
        end

        function T=ColorTable(cellTypes, sizes)
            rg=[min(sizes) max(sizes)];
            N=length(cellTypes);
            clrs=zeros(N,3);
            app=BasicMap.Global;
            cbn=app.colorsByName;
            hsls={};
            for i=1:N
                cellType=cellTypes{i};
                clr=cbn.get(cellType);
                if isempty(clr)
                    clrs(i,:)=Gui.HslColor(i,N);
                    hsls{end+1}=cellType;
                else
                    clrs(i,:)=clr;
                end
            end
            if length(hsls)<N && ~isempty(hsls)
                warning('%d/%d cell types have no color: %s', ...
                    length(hsls), N, StringArray.toString(hsls));
            end
            l=cell(N, 1);
            for i=1:N
                l{i, 1}=['<html>' Html.Symbol2(clrs(i,:),sizes(i), ...
                    rg, [6 10]) '</html>'];
            end
            varArgs={cellTypes;l};
            varArgs{end+1}='VariableNames';
            varArgs{end+1}={'Cell Type', 'Color'};
            T=table(varArgs{:});
        end

        function l=ColorList(cellTypes)
            N=length(cellTypes);
            clrs=zeros(N,3);
            app=BasicMap.Global;
            cbn=app.colorsByName;
            hsls={};
            for i=1:N
                cellType=cellTypes{i};
                clr=cbn.get(cellType);
                if isempty(clr)
                    clrs(i,:)=Gui.HslColor(i,N);
                    hsls{end+1}=cellType;
                else
                    clrs(i,:)=clr;
                end
            end
            if length(hsls)<N && ~isempty(hsls)
                warning('%d/%d cell types have no color: %s', ...
                    length(hsls), N, StringArray.toString(hsls));
            end
            l=cell(N, 1);
            for i=1:N
                clr=clrs(i,:);
                l{i}=['<font  ' Gui.HtmlHexColor(clr)...
                '>&bull;</font>'];            
            end
        end

        function [T, dataSet]=ComparisonsTable(resultsT, dataSet, ...
                testCase, f1Classifier, colorTable)
            COLS=ClassificationTable.COLUMNS;
            %backward compatibility before March 1, 2023
            hasSimilarities=size(resultsT, 2)>=8;
            col2=resultsT{:, COLS{2}};
            if isscalar(dataSet) % MLP paper clue
                isTensorFlowDataset=strcmpi(col2, 'GHOSN') ...
                    | strcmpi(col2, 'OMIP-058') ...
                    | strcmpi(col2, 'LEIPOLD');
                if strcmp(dataSet, '4')
                    dataSet='EPP''s top 4 datasets';
                    isFitcnetDataset=strcmpi(col2, 'OMIP-044') ...
                        | strcmpi(col2, 'OMIP-047') ...
                        | strcmpi(col2, 'OMIP-077') ...
                        | strcmpi(col2, 'GENENTECH');
                    isTensorFlowDataset=false(length(isFitcnetDataset), 1);
                elseif strcmp(dataSet, '6')
                    dataSet='6 datasets';
                    isFitcnetDataset=strcmpi(col2, 'OMIP-044') ...
                        | strcmpi(col2, 'OMIP-047') ...
                        | strcmpi(col2, 'GENENTECH');
                elseif strcmp(dataSet, '*')
                    dataSet='9 datasets';
                    isFitcnetDataset=strcmpi(col2, 'OMIP-044') ...
                        | strcmpi(col2, 'OMIP-047') ...
                        | strcmpi(col2, 'OMIP-069') ...
                        | strcmpi(col2, 'OMIP-077') ...
                        | strcmpi(col2, 'GENENTECH') ...
                        | strcmp(col2, 'PANORAMA');                    
                else
                    T=[];
                    return;
                end
                isNotBackground=~strcmpi(resultsT{:, COLS{4}}, ...
                    'Background');
                isTestCase=strcmpi(resultsT{:, COLS{3}}, testCase);
                col1=resultsT{:,COLS{1}};
                isFitcnet=strcmpi(col1, 'fitcnet');
                isTensorFlow=strcmpi(col1, 'TensorFlow');
                isNotMlp=~(isFitcnet|isTensorFlow);
                isDataset=isFitcnetDataset | isTensorFlowDataset ;
                isNotMlpDataset=isDataset & isNotMlp;
                isMlpDataset=(isFitcnet&isFitcnetDataset) | ...
                    (isTensorFlow&isTensorFlowDataset);
                copyT=resultsT;%safety for future programming
                copyT{isMlpDataset, COLS{1}}={'MLP'};
                T2=copyT((isNotMlpDataset | isMlpDataset)...
                    & isNotBackground...
                    & isTestCase, :);
                ds2=unique(T2{:,COLS{2}});
                cts=zeros(1, length(ds2));
                for ii=1:length(ds2)
                    ct=unique(T2{strcmp(T2{:, COLS{2}}, ds2(ii)), COLS{4}});
                    cts(ii)=length(ct);
                end
                dataSet=[num2str(sum(cts)) ' populations in ' dataSet];
            else
                T2=resultsT(strcmpi(col2, dataSet) & strcmpi(resultsT{:, COLS{3}}, testCase), :);
                ds2=unique(T2{:,COLS{2}});
                cts=zeros(1, length(ds2));
                for ii=1:length(ds2)
                    ct=unique(T2{strcmp(T2{:, COLS{2}}, ds2(ii)), COLS{4}});
                    if any(contains(ct, 'Background'))
                        cts(ii)=length(ct)-1;
                    else
                        cts(ii)=length(ct);
                    end
                end
                dataSet=[num2str(sum(cts)) ' populations in ' dataSet];
            end
            rowSize=size(T2,1);
            if rowSize<1
                T=[];
                u=unique(col2);
                [~,I]=sort(upper(u));
                html=Html.ToList(u(I));
                msgWarning(sprintf(['<html>Dataset "<font color="red">' ...
                    '%s</font>" not found!' ...
                    '<br>The datasets in this file are %s<hr></html>'], dataSet, html));
                return;
            end
            uClassifiers=StringArray.Sort(unique(T2{:,COLS{1}}), ...
                {'MLP', 'fitcnet', 'TensorFlow', 'LDA', ...
                'PhenoGraph', 'FlowSOM',  'EPP'});

            if strcmp(f1Classifier, '*')
                f1Classifier=uClassifiers{1};
            elseif ~contains(uClassifiers, f1Classifier)
                error('Classifier "%s" not found!', f1Classifier);
            end
            nClassifiers=length(uClassifiers);
            classifiers={f1Classifier};
            for i=1:nClassifiers
                if ~strcmp(uClassifiers{i}, f1Classifier)
                    classifiers{end+1}=uClassifiers{i};
                end
            end
            T3=T2(strcmp(T2{:, COLS{1}}, f1Classifier), {COLS{4}, COLS{5}, COLS{6}});
            [~, f1Order]=sort(T3{:, COLS{6}}, 'descend');
            cellTypes=T3{f1Order, COLS{4}};
            szs=T3{f1Order, COLS{5}};
            uu=unique(T2{:,COLS{4}});
            fit=ismember(uu, cellTypes);
            if ~all(fit)
                misFit=uu(~fit);
                nMisFit=length(misFit);
                for i=1:nMisFit
                    mf=misFit{i};
                    sz=T2{strcmp(mf, T2{:,COLS{4}}), COLS{5}};
                    szs(end+1)=sz(1);
                    cellTypes{end+1}=mf;
                end
            end
            l=strcmpi(cellTypes, 'background');
            if any(l)
                cellTypes(l)=[];
                szs(l)=[];
            end
            nCellTypes=length(cellTypes);
            if nargin>4
                if ~isa(colorTable, 'table')
                    colorTable=ClassificationTable.ColorTable(cellTypes, ...
                        double(szs));
                end
                bullets=cell(nCellTypes, 1);
                for i=1:nCellTypes
                    cellType=cellTypes{i};
                    c=colorTable{strcmp(colorTable{:, 1}, cellType), 2};
                    if ~isempty(c)
                        bullets{i}=c{1};
                    else
                        bullets{i}=['<html><font size="6"><font ' ...
                            'color="#777777">&bull;</font></font>' ...
                            '</</html>'];
                    end
                end
                varArgs={bullets;cellTypes;szs};
                cols={'Color', 'Cell type', 'Size'};
            else
                cols={'Cell type', 'Size'};
                varArgs={cellTypes;szs};
            end
            f1s=nan(nCellTypes, nClassifiers);
            ics=nan(nCellTypes, nClassifiers);
            similarities=nan(nCellTypes, nClassifiers);
            for i=1:nCellTypes
                cellType=cellTypes{i};
                l1=strcmp(T2{:, COLS{4}}, cellType);
                for j=1:nClassifiers
                    classifier=classifiers{j};
                    l2=strcmp(T2{:, COLS{1}}, classifier);
                    if any(l1 & l2)
                        f1=T2{l1 & l2, COLS{6}};
                        ic=T2{l1 & l2, COLS{7}};
                        f1s(i,j)=f1(1);
                        ics(i,j)=ic(1);
                        if hasSimilarities
                            similarity=T2{l1 & l2, COLS{8}};
                            similarities(i,j)=similarity(1);
                        end
                    end
                end
            end
            for j=1:nClassifiers
                cols{end+1}=[classifiers{j} ' F1'];
                varArgs{end+1}=f1s(:, j);
            end
            for j=1:nClassifiers
                cols{end+1}=[classifiers{j} ' IC'];
                varArgs{end+1}=ics(:, j);
            end
            if hasSimilarities
                for j=1:nClassifiers
                    cols{end+1}=[classifiers{j} ' Similarity'];
                    varArgs{end+1}=similarities(:, j);
                end
            end
            varArgs{end+1}='VariableNames';
            varArgs{end+1}=cols;
            T=table(varArgs{:});
        end
    
        function found=Find(T, cellTypes, classifier, dataset, testCase)
            N=length(cellTypes);
            found=zeros(N,1);
            COLS=ClassificationTable.COLUMNS;
            l=TableBasics.Select(T, COLS{1}, classifier, [], '=', 'rows');
            if any(l)
                l=l&TableBasics.Select(T, COLS{2}, dataset, [], '=', 'rows');
                if any(l)
                    l=l&TableBasics.Select(T, COLS{3}, testCase, [], '=', 'rows');
                    if any(l)
                        if isnumeric(cellTypes)
                            nCellTypes=length(cellTypes);
                            strs=cell(1, nCellTypes);
                            for i=1:nCellTypes
                                strs{i}=num2str(cellTypes(i));
                            end
                            cellTypes=strs;
                        end
                        a=T{:, 'CellType'};
                        for i=1:N
                           l0=strcmpi(a, cellTypes{i});
                           idx=find(l&l0, 1);
                           if ~isempty(idx)
                               found(i)=idx;
                           end
                        end
                    end
                end
            end
        end

        function [T, file]=Get(ask, jw, propFile, propFldr, props, folder)
            if nargin<6
                folder=File.Documents;
                if nargin<5
                    props=BasicMap.Global;
                    if nargin<4
                        propFldr=ClassificationTable.PROP_FOLDER;
                        if nargin<3
                            propFile=ClassificationTable.PROP_FILE;
                            if nargin<2
                                jw=[];
                                if nargin<1
                                    ask=false;
                                end
                            end
                        end
                    end
                end
            end
            if isa(jw, 'QfTable')
                jw=Gui.JWindow(jw.fig);
            elseif isstruct(jw) && isfield(jw, 'jw')
                jw=jw.jw;
            else
                jw=Gui.JWindow(get(0, 'CurrentFigure'));
            end
            T=[];
            file=[];
            while true
                file=TableBasics.UiPutFile(propFldr, propFile, ...
                    props, '5_classifiers.xls', folder);
                if isempty(file) || ~exist(file, 'file')
                    return;
                else
                    if ~endsWith(file, 'xls')
                        T=readtable(file, 'Delimiter', ',');
                    else
                        T=readtable(file);
                    end
                end
                COLS=ClassificationTable.COLUMNS;
                l=ismember(COLS, T.Properties.VariableNames);
                if ~all(l)
                    COLS=COLS(~l);
                    if ~ask
                        error('File is missing columns %s', ...
                            StringArray.toString(COLS, ',',true));
                    else
                        html=Html.ToList(COLS);
                        if ~askYesOrNo(struct('javaWindow', jw, 'msg', ...
                                ['<html>File is missing columns' html ...
                                '<b>Try again?</b><hr></html>']))
                            T=[];
                            file=[];
                            return;
                        end
                    end
                else
                    return;
                end
            end
        end

        function [T, file, classifierName, datasetName, testCaseName, results, viewToo]...
                =IntegrateResults(results, ask, T, file, classifierName,...
                datasetName, testCaseName, doSortQF)
            if verLessThan('matlab', '9.11')
                msgError(['<html><center>MATLAB version R2021b '...
                '<br>or later is needed for this</center></html>']);
                T=[];
                file=[]; classifierName=[]; datasetName=[];
                testCaseName=[]; results=[];
                return;
            end
            if nargin<8
                doSortQF=true;
                if nargin<7
                    testCaseName=[];
                    if nargin<6
                        datasetName='';
                        if nargin<5
                            classifierName='';
                            if nargin<4
                                file='';
                                if nargin<3
                                    T=[];
                                    if nargin<2
                                        ask=false;
                                        if nargin<1
                                            error(['Need a source with ' ...
                                                'which to integrate...']);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if isempty(testCaseName)
                testCaseName='all samples';
            end
            isQft=isa(results, 'QfTable');
            if isQft
                qfT=results.sortTable.getSortedTableData(results.data);
                cols={'Subset (class) name', ...
                    '# of events', ...
                    '% Overlap (F1-Score)', ...
                    'Central similarity (inlier concordance)', ...
                    'Similarity, earth mover''s distance', ...
                    '% Concordance (Jaccard index)'};
                qfT=TableBasics.Select(qfT, 'Training set?', 'yes', cols,[], 'table');
                results=struct();
                if ~doSortQF
                    results.cellTypes=qfT{:, 1};
                    results.size=qfT{:, 2};
                    results.f1=qfT{:, 3};%/100;
                    results.ic=qfT{:, 4};%/100;
                    results.similarity=qfT{:, 5};
                    results.ji=qfT{:, 6};
                else
                    ST=qfT{:, 1};
                    [~, I]=sort(upper(ST));
                    results.cellTypes=ST(I);
                    results.size=qfT{I, 2};
                    results.f1=qfT{I, 3};%/100;
                    results.ic=qfT{I, 4};%/100;
                    results.similarity=qfT{I, 5};%/100;
                    results.ji=qfT{I, 6};
                end
            end
            viewToo=true;
            if ask
                [cancelled, file, datasetName, classifierName, viewToo]=...
                    ClassificationTable.Ask(datasetName, ...
                    classifierName, file);
                if cancelled
                    T=[];
                    file=[];
                    return;
                end
            end
            if ~isempty(file)
                try
                    T=readtable(file);
                catch
                end
            elseif isempty(T)
                [T, file]=ClassificationTable.Get(ask, results);
                if isempty(file)
                    T=[];
                    file=[];
                    return;
                end
            end
            R=length(results.cellTypes);
            if isempty(T)
                T=[];
                needToAdd=true(1,R);
                needToUpdate=false;
            else
                found=ClassificationTable.Find(T, results.cellTypes, ...
                    classifierName, datasetName, testCaseName);
                needToAdd=found==0;
                needToUpdate=~needToAdd;
                updateIdxs=found(needToUpdate);
            end
            if size(T, 2)==7
                Similarity=zeros(size(T,1),1);
                TT=table(Similarity);
                T=[T TT];
            end
            COLS=ClassificationTable.COLUMNS;
            R=sum(needToAdd);
            if R>0
                newT=ClassificationTable.New(R);
                newT{:, COLS{1}}={classifierName};
                newT{:, COLS{2}}={datasetName};
                newT{:, COLS{3}}={testCaseName};
                newT{:, COLS{4}}=results.cellTypes(needToAdd);
                newT{:, COLS{5}}=results.size(needToAdd);
                if isfield(results, 'f1')
                    newT{:, COLS{6}}=results.f1(needToAdd);
                end
                if isfield(results, 'ic')
                    newT{:, COLS{7}}=results.ic(needToAdd);
                end
                if isfield(results, 'similarity')
                    newT{:, COLS{8}}=results.similarity(needToAdd);
                end
                if isfield(results, 'ji')
                    newT{:, COLS{9}}=results.ji(needToAdd);
                end

                T=[T;newT];
            end
            if any(needToUpdate)
                T{updateIdxs, COLS{1}}={classifierName};
                T{updateIdxs, COLS{2}}={datasetName};
                T{updateIdxs, COLS{3}}={testCaseName};
                T{updateIdxs, COLS{4}}=results.cellTypes(needToUpdate);
                T{updateIdxs, COLS{5}}=results.size(needToUpdate);
                if isfield(results, 'f1')
                    T{updateIdxs, COLS{6}}=results.f1(needToUpdate);
                end
                if isfield(results, 'ic')
                    T{updateIdxs, COLS{7}}=results.ic(needToUpdate);
                end
                if isfield(results, 'similarity')
                    T{updateIdxs, COLS{8}}=results.similarity(needToUpdate);
                end
                if isfield(results, 'ji')
                    T{updateIdxs, COLS{9}}=results.ji(needToUpdate);
                end
            end
            if ~isempty(file)
                writetable(T, file);
            end
        end

        function [T, figBrowser, f1Fig, icFig]=...
                IntegrateLdaMlp(f1T, icT, simT, jiT, ...
                dataSet, show, file, ask, mlpName)
            if nargin<9
                mlpName='fitcnet';
                if nargin<8
                    ask=true;
                    if nargin<7
                        file=[];
                        if nargin<6
                            show=2;
                        end
                    end
                end
            end
            figBrowser=[];f1Fig=[]; icFig=[];
            testCase='all samples';
            resultsLda=struct();
            cellTypes=f1T{:, 1};
            if isnumeric(cellTypes)
                nCellTypes=length(cellTypes);
                strs=cell(nCellTypes, 1);
                for i=1:nCellTypes
                    strs{i, 1}=num2str(cellTypes(i));
                end
                cellTypes=strs;
            end
            resultsLda.cellTypes=cellTypes;
            resultsLda.f1=f1T{:, 2};
            resultsLda.ic=icT{:, 2};
            resultsLda.similarity=simT{:,2};
            resultsLda.ji=jiT{:,2};
            resultsLda.size=f1T{:, 6};
            resultsMlp=struct();
            resultsMlp.cellTypes=cellTypes;
            resultsMlp.f1=f1T{:, 3};
            resultsMlp.ic=icT{:, 3};
            resultsMlp.similarity=simT{:,3};
            resultsMlp.ji=jiT{:,3};
            resultsMlp.size=f1T{:, 7};
            if ~isempty(file)
                if ~exist(file, 'file')
                    file2=ClassificationTable.DownloadClasses;
                    if ~isempty(file2)
                        copyfile(file2, file);
                    end
                end
            end
            [T, tFile, ~, dataSet]=ClassificationTable.IntegrateResults( ...
                resultsMlp, ask, [], file, mlpName, dataSet, testCase);
            if ~isempty(T)
                T=ClassificationTable.IntegrateResults(resultsLda, ...
                    false, T, tFile, 'LDA', dataSet, testCase);
            end
            if ~isempty(T)
                [~,~, figBrowser,f1Fig,icFig]=ClassificationTable.Html(...
                    ClassificationTable.ComparisonsTable( ...
                    T, dataSet, testCase, mlpName, true), ...
                    dataSet, show);
            end
        end

        function [whichMetrics, cancelled]=ChoosePerformanceMetrics
            [whichMetrics, cancelled]=Gui.Ask('Which performance metrics?', {...
                'F1-score and CS (AKA JI 80%)', ...
                'JI 100% and JI 80%', ...
                'F1-Score, CS and EMD', ...
                'F1-score, JI 100%, JI-80% and EMD',...
                'F1-score ONLY'}, ...
                'MlpPaperFunctions.WhichMetric', 'Performance metric');
            if cancelled
                return;
            end
            if whichMetrics==3
                whichMetrics=0;
            end
        end

        function [f1Fig, icFig, emdFig]=Plots(comparisonTable, ...
                dataSetName, show, followed, where, plotMap, ...
                hasSimilarities, args)
            if nargin<8
                args=struct(...
                    'plotMdnMn', true, ...
                    'inlierConcordanceName', false, ...
                    'similarityName', false, ...
                    'genentechName', 'ESHGHI');
            end
            if args.inlierConcordanceName
                icName='Inlier Concordances';
            else
                icName='Central similarity';
            end
            if args.similarityName
                simName='Similarity';
            else
                simName='Earth mover''s distance';
            end
            [f1Fig, tb]=Gui.NewFigure(true, false, true);
            set(f1Fig, 'Name', ['F1-scores for ' dataSetName]);
            set(f1Fig, 'Color', 'white')
            Gui.AddSvgToToolBar(f1Fig, tb);
            [icFig, tb]=Gui.NewFigure(true, false, true);
            set(icFig, 'Color', 'white')
            set(icFig, 'Name', [icName ' for ' dataSetName])
            Gui.AddSvgToToolBar(icFig, tb);
            C=size(comparisonTable ,2);
            if hasSimilarities
                [emdFig, tb]=Gui.NewFigure(true, false, true);
                set(emdFig, 'Color', 'white')
                set(emdFig, 'Name', [simName ' for ' dataSetName])
                Gui.AddSvgToToolBar(emdFig, tb);
            else
                emdFig=[];
            end
            if isempty(comparisonTable.Properties.VariableDescriptions)
                H=comparisonTable.Properties.VariableNames;
            else
                H=comparisonTable.Properties.VariableDescriptions;
            end
            classifiers={};
            clrs=comparisonTable{:, 1};
            CellTypes=comparisonTable{:, 2};
            sizes=comparisonTable{:, 3};
            ics=[];
            f1s=[];
            similarities=[];
            for c=1:C
                h=H{c};
                if endsWith(h, ' IC')
                    H{c}=h(1:end-3);
                    ics=[ics comparisonTable{:, h}];
                elseif endsWith(h, ' F1')
                    H{c}=h(1:end-3);
                    classifiers{end+1}=H{c};
                    f1s=[f1s comparisonTable{:, h}];
                elseif endsWith(h, ' Similarity')
                    H{c}=h(1:end-11);
                    similarities=[similarities comparisonTable{:, h}];
                end
            end
            nClassifiers=length(classifiers);
            [plotRows, plotCols]=Gui.GetSubPlotSize(nClassifiers);
            for c=1:nClassifiers
                classifier=classifiers{c};
                f1Ax=subplot(plotRows, plotCols, c, 'Parent', f1Fig);
                axis(f1Ax, 'square');
                scores=f1s(:, c);
                [~, ~, ~, ~, clrs]=Plots.ScoreSizes([], scores, sizes,...
                    CellTypes, classifier, 'F1-Score', [], [], ...
                    clrs, f1Ax, plotMap, args);
            end
            for c=1:nClassifiers
                classifier=classifiers{c};
                icAx=subplot(plotRows, plotCols, c, 'Parent', icFig);
                axis(icAx, 'square');
                scores=ics(:, c);
                Plots.ScoreSizes([], scores, sizes, CellTypes, ...
                    classifier, 'Inlier Concordance', [], [], ...
                    clrs, icAx, plotMap, args);
            end
            if hasSimilarities
                for c=1:nClassifiers
                    classifier=classifiers{c};
                    simAx=subplot(plotRows, plotCols, c, 'Parent', emdFig);
                    axis(simAx, 'square');
                    scores=similarities(:, c);
                    Plots.ScoreSizes([], scores, sizes, CellTypes, ...
                        classifier, 'Similarity', [], [], ...
                        clrs, simAx, plotMap, args);
                end
            end
            if nargin>3
                SuhWindow.Follow(f1Fig, followed, where, true);
            end
            app=BasicMap.Global;
            finishFig(f1Fig, 'F1-score', show);
            SuhWindow.Follow(icFig, f1Fig, 'south++', true);
            finishFig(icFig, icName, show);
            if hasSimilarities
                SuhWindow.Follow(emdFig, icFig, 'west++', true);
                finishFig(emdFig, simName, show);
            end

            function finishFig(fig, scoreName, show)
                han=axes(fig,'visible','off');
                axis(han, 'square');
                set(get(han, 'Title'), 'Visible','on');
                set(get(han, 'XLabel'), 'Visible', 'on');
                set(get(han, 'YLabel'), 'Visible', 'on');
                yl=ylabel(han, scoreName);
                xl=xlabel(han,'Number of cells in population (log10)');
                clr=ClassificationTable.COLOR;
                if app.highDef
                    nudge=1;
                else
                    nudge=3;
                end
                set(yl, 'FontSize', get(yl, 'FontSize')+nudge, 'Color', clr);
                set(xl, 'FontSize', get(xl, 'FontSize')+nudge, 'Color', clr);
                try
                    if app.highDef
                        fs=13;
                    else
                        fs=14;
                    end
                    if args.plotMdnMn
                        sgtitle(fig, [fig.Name ' median/mean/weighted'], ...
                            'FontSize', fs, 'Color', clr);
                    else
                        sgtitle(fig, fig.Name , ...
                            'FontSize', fs, 'Color', clr);
                    end
                    set(fig, 'Name', String.RemoveTex(fig.Name));
                catch ex
                    ex.getReport
                end
                if show
                    SuhWindow.SetFigVisible(fig);
                end
                yp=get(yl, 'Position');
                if plotRows<plotCols
                    yp(1)=yp(1)-(plotRows/10);
                    if app.highDef
                        yp(1)=yp(1)*.8;
                    end
                else
                    yp(1)=yp(1)*2;
                end
                set(yl, 'Position', yp);
            end
        end


        function [html, header, figBrowser, f1Fig, icFig, emdFig]=Html( ...
                comparisonTable, dataSetName, ...
                showPlots, followed, where, fileToEdit, args)
            if nargin<7
                args=struct(...
                    'plotMdnMn', false, ...
                    'inlierConcordanceName', false, ...
                    'similarityName', false, ...
                    'genentechName', 'ESHGHI');                
            end
            if strcmpi('genentech', dataSetName)
                dataSetName=args.genentechName;
            end
            if isempty(comparisonTable)
                msgWarning('No comparisons found ....');
                html=[];header=[]; figBrowser=[]; 
                f1Fig=[]; icFig=[]; emdFig=[];
                return;
            end
            app=BasicMap.Global;
            f1Fig=[]; icFig=[];emdFig=[];
            if nargin<6
                fileToEdit=[];
                if nargin<5
                    where='east+';
                    if nargin<4
                        followed=[];
                        if nargin<3
                            showPlots=true;
                        end
                    end
                end
            end
            totalSize=sum(comparisonTable{:,3});
            [R,C]=size(comparisonTable);
            %MUST remove xml from cell type name
            CellTypes=comparisonTable{:, 2};
            sortableCellTypes=cell(1, R);
            for r=1:R
                v=Html.Remove(CellTypes{r});
                sortableCellTypes(1,r)=upper( ...
                    edu.stanford.facs.swing.MarkerSorter.encodeKey(v));
                CellTypes{r}=v;
            end
            try
                comparisonTable(:,2)=CellTypes;
            catch
                comparisonTable(:,2)=table(CellTypes);
            end
            if isempty(comparisonTable.Properties.VariableDescriptions)
                header=comparisonTable.Properties.VariableNames;
            else
                header=comparisonTable.Properties.VariableDescriptions;
            end
            hasSimilarities=false;
            classifiers={};
            for c=1:C
                h=header{c};
                if endsWith(h, ' IC')
                    header{c}=h(1:end-3);
                elseif endsWith(h, ' F1')
                    header{c}=h(1:end-3);
                    classifiers{end+1}=header{c};
                elseif endsWith(h, ' Similarity')
                    hasSimilarities=true;
                    header{c}=h(1:end-11);
                end
            end
            nClassifiers=length(classifiers);
            if nClassifiers==1
                where2='west';
                if ~isempty(followed)
                    if strcmpi(where, 'west')
                        where2='north';
                    end
                end
                msg(Html.SprintfHr(['No other classifications found ' ...
                    'for<br>dataset "<b>%s</b>."'], dataSetName), 6, where2);
            end
            sortCol=4;
            if hasSimilarities
                ascending=false(3+nClassifiers+nClassifiers+nClassifiers,1);
            else
                ascending=false(3+nClassifiers+nClassifiers,1);
            end
            ascending(4)=true;%f1 score
            theTitle=sprintf(['%s scores for %s '],...
                Html.Remove(String.RemoveTex(dataSetName)), ...
                String.Pluralize2('classifier', nClassifiers));
            html=doHtml;
            [figBrowser, tp, btnForTip, ~]=internalBrowser(html);
            %shift window down in size a tad to position as legend
            pos=Gui.GetOuterPixels(figBrowser);
            dim=tp.getPreferredSize;
            if dim.width<pos(3)
                pos(3)=dim.width;
            else
                pos(3)=pos(3)*.80;
            end
            if dim.height<pos(4)
                pos(4)=dim.height*1.02;
            else
                pos(4)=pos(4)*.75;
            end
            set(figBrowser, 'OuterPosition', pos);
            SuhWindow.SetFigVisible(figBrowser);
            if ~isempty(followed)
                SuhWindow.Follow(figBrowser, followed, where);
            end
            if showPlots
                plotMap=Map;
                keys={};
                for i=1:nClassifiers
                    keys{end+1}=[classifiers{i} '.F1-Score.'];
                    keys{end+1}=[classifiers{i} '.Inlier Concordance.'];
                    if hasSimilarities
                        keys{end+1}=[classifiers{i} '.Similarity.'];
                    end
                end
                nKeys=length(keys);
                [f1Fig, icFig, emdFig]=ClassificationTable.Plots( ...
                    comparisonTable, dataSetName, true, figBrowser, ...
                    'north east++', plotMap, hasSimilarities, args);
            end
            
            function html=doHtml
                sb=java.lang.StringBuilder();
                sb.append("<br><table border='2' cellpadding='0' cellspacing='0'><thead><tr>");
                for col=1:3
                    sb.append("<th rowspan='2'>&nbsp;");
                    if col==1
                        sb.append(header{col});
                    else
                        if col==3
                            doSortHeader(sb, col, ['Frequ<br>' ...
                                '&nbsp;&nbsp;&nbsp;ency %']);
                        else
                            doSortHeader(sb, col, header{col});
                        end
                    end
                    sb.append("&nbsp;</th>");
                end
                if nClassifiers<3
                    br='<br>';
                else
                    br=' ';
                end
                if args.inlierConcordanceName
                    inlierConcordanceName=['<th colspan=' num2str(nClassifiers)...
                        '" align="center">&nbsp;Inlier' br ...
                        'concordance&nbsp;</th>'];
                else
                    inlierConcordanceName=['<th colspan=' num2str(nClassifiers)...
                        '" align="center">&nbsp;Central' br ...
                        'similarity&nbsp;</th>'];
                end
                if hasSimilarities
                    if args.similarityName
                        sb.append(['<th colspan="' num2str(nClassifiers)...
                            '" align="center">&nbsp;F1-Score&nbsp;</th>' ...
                            inlierConcordanceName '<th colspan="' ...
                            num2str(nClassifiers)...
                            '" align="center">&nbsp;Similarity&nbsp;'...
                            '</th></tr><tr>']);
                    else
                        sb.append(['<th colspan="' num2str(nClassifiers)...
                            '" align="center">&nbsp;F1-Score&nbsp;</th>'...
                            inlierConcordanceName '<th colspan="'...
                            num2str(nClassifiers)...
                            '" align="center">&nbsp;Earth mover''s' ...
                            ' distance&nbsp;</th></tr><tr>']);
                    end
                else
                    sb.append(['<th colspan="' num2str(nClassifiers)...
                        '" align="center">&nbsp;F1-Score&nbsp;</th>' ...
                        inlierConcordanceName '</tr><tr>']);
                end
                classHeaders(sb, 4);
                classHeaders(sb, 4+nClassifiers);
                if hasSimilarities
                    classHeaders(sb, 4+(2*nClassifiers));
                end
                sb.append(sprintf("</tr></thead>"));
                if R>0
                    conv=cell(1,C);
                    for col=1:C
                        value=comparisonTable{1,col};
                        if isnumeric(value)
                            conv{col}=@nCell;
                        else
                            conv{col}=@sCell;
                        end
                    end
                    for r=1:R
                        sb.append("<tr>");
                        for col=1:C
                            value=comparisonTable{r,col};
                            feval(conv{col}, value, col, sb);
                        end
                        sb.append(sprintf("</tr>"));
                    end
                end
                try
                    sb.append(TableBasics.Summarize(comparisonTable, ...
                        {'Median', 'Mean'}, {@median, @mean}, 3, ...
                        4:(3+(nClassifiers*3)), @nCellSummarize,...
                        {"<tr><th colspan='2' rowspan='2'>" + num2str(R) + ...
                        " cell populations</th>", ""}, 1));
                catch ex
                    ex.getReport
                end
                sb.append("</table>");
                html=char(sb.toString);
            end

            function classHeaders(sb, col)
                for cc=1:nClassifiers
                    sb.append("<th align='center'>&nbsp;");
                    classifier=classifiers{cc};
                    nCl=length(classifier);
                    if nCl>6
                        half=ceil(nCl/2);
                        classifier=[classifier(1:half) ...
                            '<br>&nbsp;' classifier(half+1:end)];
                    end
                    doSortHeader(sb, col+cc-1, classifier);
                    sb.append("</th>");
                end
            end

            function doSortHeader(sb, theCol, value)
                sb.append(sprintf("<a href='%d'>", theCol));
                if theCol==sortCol
                    if ~ascending(sortCol)
                        sb.append("<font color='green'>");
                    else
                        sb.append("<font color='#229999'>");
                    end
                    sb.append('<i>');
                    sb.append(value);
                    sb.append('</i>');
                    sb.append("</a>&nbsp;");
                    sb.append("</font>");
                    if ~ascending(sortCol)
                        sb.append(Html.Img('sortDown.png'));
                    else
                        sb.append(Html.Img('sortUp.png'));
                    end
                else
                    sb.append(value);
                    sb.append("</a>&nbsp;");
                end

            end

            function out=browserHtml
                ss=strrep(html, "<tr>", sprintf("\n<tr>"));
                ss=strrep(ss, 'file:/%2F', 'file:/');
                ss=strrep(ss, '%2F', '/');
                out=Html.Wrap(...
                    ['<h2>' String.RemoveTex(theTitle) '</h2>' ...
                    char(ss)]);
            end

            function flash(cellType)
                figure(f1Fig);
                figure(icFig);
                if ~isempty(emdFig)
                    figure(emdFig);
                end
                drawnow;
                szs=zeros(1, nKeys);
                plotHs=zeros(1, nKeys);
                for k=1:nKeys
                    try
                        key=[keys{k} cellType];
                        plotHs(k)=plotMap.get(key);
                        szs(k)=get(plotHs(k), 'MarkerSize');
                    catch ex
                        warning('No plot for %s?\n\t%s', key, ex.message);
                    end
                end
                for kk=1:3
                    for k=1:nKeys
                        if szs(k)>0
                            set(plotHs(k), 'MarkerSize', szs(k)*2);
                            nSecs=Gui.FlashN(plotHs(k), 1, .06, false, true);
                        end
                    end
                end
                MatBasics.RunLater(@(h,e)resizes(plotHs, szs), nSecs*2);
            end

            function resizes(plotHs, szs)
               for k=1:nKeys
                   resize(plotHs(k), szs(k));
               end
            end

            function resize(H, sz)
                set(H, 'MarkerSize', sz);
            end

            function [figBrowser, tp, btnForTip, btnSvg]=internalBrowser(html)
                [figBrowser, tp, ~, btnSvg, tb]=Gui.FigBrowser( ...
                    html, @hyper, @browserHtml, false);
                set(figBrowser, 'Name', theTitle);
                if ~isempty(fileToEdit)
                    btnForTip=ToolBarMethods.addButton(tb, 'table.gif', ...
                        ['Open the classification file ' ...
                        'supporting this view'], ...
                        @(h,e)openClassification(),...
                        'See underlying data file');
                else
                    btnForTip=btnSvg;
                end
            end

            function openClassification
                [choice, cancelled]=Gui.Ask(['<html>' ...
                    Html.FileTree(fileToEdit) ...
                    '<br><br>How do wish to open this file?</html>'], ...
                    {'Open the classification file DIRECTLY', ...
                    'Open the folder for the classification file'}, ...
                    'ClassificationTable.Html.Open', ...
                    'See underlying data....', 1);
                if cancelled
                    return;
                end
                if choice==1
                    if ismac
                        system(['open ' String.ToSystem(fileToEdit)]);
                    else
                        system(String.ToSystem(fileToEdit));
                    end
                else
                    File.OpenFolderWindow(fileToEdit,[], false);
                end
                MatBasics.RunLater(@(h,e)warnToSave, 2);
            end

            function warnToSave
                msgWarning(struct('javaWindow', Gui.JWindow(figBrowser),...
                    'msg', Html.WrapHr(['If you save changes you must' ...
                    '<br>re-generate the table and plots.' ...
                    '<br><br>' Html.WrapSmallBoldOnly(['(Also: ' ...
                    'close quickly if using <i>Microsoft</i> Excel)'])])), 12, 'center', ...
                    'Sorry but ...');
            end

            function hyper(~, eventData)
                tipX=0; tipY=30;
                description = char(eventData.getDescription); % URL stri
                et=char(eventData.getEventType);
                switch char(et)
                    case char(eventData.getEventType.ENTERED)
                        col=str2num(description);
                        if isempty(col) || isnan(col)
                            app.showToolTip(btnForTip, ['<html>Flash <b>' ...
                                description '</b> in the plots...'], ...
                                tipX, tipY);
                        else
                            if col>3
                                idx=col-3;
                                if idx>nClassifiers*2
                                    idx=idx-(nClassifiers*2);
                                    score='similarity';
                                elseif idx>nClassifiers
                                    idx=idx-nClassifiers;
                                    score='inlier concordance';
                                else
                                    score='F1-score';
                                end
                                classifier=classifiers{idx};
                            else
                                classifier='this column';
                                if col==3
                                    score='size';
                                elseif col==2
                                    score='name value';
                                else
                                    score='';
                                end
                            end
                            if ascending(col)
                                app.showToolTip(btnForTip, ['<html>Sort ' ...
                                    'table by <i>ascending</i> ' ...
                                    score ' for <b>' classifier ...
                                    '</b></html>'], ...
                                    tipX, tipY);
                            else
                                app.showToolTip(btnForTip, ['<html>Sort ' ...
                                    'table by <i>descending</i> ' ...
                                    score ' for <b>' classifier ...
                                    '</b></html>'], ...
                                    tipX, tipY);
                            end
                        end
                    case char(eventData.getEventType.EXITED)
                        %disp('link hover exit');
                    case char(eventData.getEventType.ACTIVATED)
                        fprintf('Activated "%s"\n', description);
                        if isempty(plotMap)
                            disp('NO PLOTS ....');
                            return;
                        end
                        col=str2num(description);
                        if isempty(col) || isnan(col)
                            if startsWith(description, 'cell:')
                                flash(description(6:end));
                            else
                                flash(description);
                            end
                        else
                            Gui.ShowBusy(figBrowser, Gui.YellowH2(...
                                'Sorting table'), 'genie.png', 1.1)
                            values=comparisonTable{:, col};
                            if col ~= 2
                                values(isnan(values))=0;
                                if ascending(col)
                                    [~, order]=sort(values, 'ascend');
                                else
                                    [~, order]=sort(values, 'descend');
                                end
                            else
                                % case insensitive sort on name
                                [~, order]=...
                                    sort(sortableCellTypes);
                                if ~ascending(col)
                                    order=flip(order);
                                end
                            end
                            sortableCellTypes=sortableCellTypes(order);
                            ascending(col)=~ascending(col);
                            comparisonTable=comparisonTable(order,:);
                            sortCol=col;
                            html=doHtml;
                            if app.highDef
                                tp.setText(['<html><body ' ...
                                    'bgcolor="#FFFFFF" style=' ...
                                    '"font-size:110%">' html ...
                                    '</body></html>']);
                            else
                                tp.setText(['<html><body bgcolor=' ...
                                    '"#FFFFFF">' html '</body></html>']);
                            end
                            tp.setCaretPosition(0)
                            Gui.HideBusy(figBrowser);
                        end
                end
            end

            function nCellSummarize(x, col, sb)
                if isnan(x)
                    sb.append('<td>&nbsp;N/A&nbsp;&nbsp;</td>');
                else
                    if col==3
                        x=ceil(x);
                    end
                    
                    if col~=3 && mod(col-3,nClassifiers)==1
                        sb.append(['<td align="right" bgcolor="#F2F2F0">' ...
                            '<b>']);
                    else
                        sb.append('<td align="right">&nbsp;&nbsp;&nbsp;<b>');
                    end
                    
                    sb.append(String.encodeRounded(x, 3, true));
                    sb.append("</b>&nbsp;</td>");
                end
            end

            function nCell(x, col, sb)
                if isnan(x)
                    sb.append('<td></td>');
                else
                    if col~=3 && mod(col-3,nClassifiers)==1
                        sb.append(['<td align="right" bgcolor="#F2F2F0">' ...
                            '&nbsp;&nbsp;&nbsp;']);
                    else
                        sb.append('<td align="right">&nbsp;&nbsp;&nbsp;');
                    end
                    if col~=3
                        if x>0
                            sb.append(String.encodeRounded(x, 2, true));
                        end
                    else
                        sb.append(String.encodePercent(double(x),totalSize));
                    end
                    sb.append("&nbsp;</td>");
                end
            end

            function sCell(x, col, sb)
                x=x{1};
                if isempty(x)
                    sb.append("<td></td>");
                elseif startsWith(x, "<html>")
                    sb.append("<td align='center'>");
                    sb.append(Html.remove(x));
                    sb.append("</td>");
                else
                    x=Html.Remove(x);
                    sb.append("<td>&nbsp;");
                    if col==2
                        sb.append('<a href="cell:');
                        sb.append(x);
                        sb.append('">');
                        sb.append(x);
                        sb.append("</a>");
                    else
                        sb.append(x);
                    end
                    sb.append("&nbsp;</td>");
                end
            end
        end

        function [cancelled, classifierName, dataSetName]=...
                GetClassifierAndDataSet(classifierName, dataSetName, ...
                prompt, where)
            if nargin<4
                where='center';
                if nargin<3
                    prompt='Add this classification to a comparison file?';
                end
            end
            cancelled=false;
            while true
                values=inputsDlg(prompt, 'Confirm...', ...
                    {'Dataset name', 'Classifier name'}, ...
                    {dataSetName, classifierName}, where, false,...
                    14, 2, 1, [],[],[],[],false,[],[],[], false);
                if isempty(values)
                    cancelled=true;
                    return;
                end
                if ~isempty(values{2})
                    classifierName=values{2};
                end
                if ~isempty(values{1})
                    dataSetName=values{1};
                end
                if ~isempty(StringArray.IndexesOfEmpties(values))
                    if ~askYesOrNo(Html.WrapHr(['No empty ' ...
                            'entries allowed ... ' ...
                            '<br>Enter again?']),[], 'center', ...
                            true,[],'ClassificationTable.Empties')
                        break
                    end
                else
                    break;
                end
            end
            if length(values)==2
                classifierName=values{2};
                dataSetName=values{1};
            else
                cancelled=true;
            end
        end

        function [T, dataSetName, figBrowser, f1Fig, icFig, emdFig]=See( ...
                dataSetName, testCaseName, file, followed, where, plotMdnMn)
            T=[];
            figBrowser=[];
            f1Fig=[];
            icFig=[];
            emdFig=[];
            if nargin<6
                plotMdnMn=true;
                if nargin<5
                    where='east+';
                    if nargin<4
                        followed=[];
                        if nargin<3
                            file=[];
                            if nargin<2
                                testCaseName=[];
                                if nargin<1
                                    dataSetName='OMIP-077';
                                end
                            end
                        end
                    end
                end
            end
            if isempty(testCaseName)
                testCaseName='all samples';
            end
            if strcmpi(dataSetName, 'ESHGHI')
                dataSetName='GENENTECH';
            end
            if isempty(file)
                [cancelled, file, dataSetName]=ClassificationTable.Ask(dataSetName);
                if cancelled
                    return;
                end
            end
            T=[];
            if isempty(file)
                file=TableBasics.UiGetFile;
            end
            pu=PopUp('Gathering classifications...', 'center', ...
                'Patience...', true, false, 'match.png', false, ...
                [], Gui.JWindow(followed) );
            if ~isempty(file)
                try
                    T=readtable(file);
                catch
                end
            end
            if isempty(T)
                pu.close;
                return;
            end
            try
                if nargin==1 && endsWith(file, '_no0s.xls')
                    T=ClassificationTable.ComparisonsTable( ...
                        T, dataSetName, testCaseName, 'MLP', true);
                else
                    [T, dataSetName]=ClassificationTable.ComparisonsTable( ...
                        T, dataSetName, testCaseName, 'fitcnet', true);
                end
            catch ex
                ex.getReport
                [T, dataSetName]=ClassificationTable.ComparisonsTable( ...
                    T, dataSetName, testCaseName, '*', true);
            end
            if ~isempty(T)
                  args=struct(...
                    'plotMdnMn', plotMdnMn, ...
                    'inlierConcordanceName', false, ...
                    'similarityName', false, ...
                    'genentechName', 'ESHGHI');      
                [~, ~, figBrowser,f1Fig,icFig, emdFig]=...
                    ClassificationTable.Html(...
                    T, dataSetName, true, followed, where, file, args);
            end
            pu.close;
        end


        function [bp, jtfDataset, jtfClassifier]=EntryPanel(dataset, classifier)
            [jtfDataset, west]=make('Dataset', dataset, ['<html>Enter the ' ...
                'dataset name<br>(e.g.OMIP-047, GHOSN, etc.)</html>']);
            if nargin<2
                jtfClassifier=[];
                bp=Gui.BorderPanel([],0,0, 'West', west);
            else
                [jtfClassifier, east]=make('Classifier', classifier,...
                    ['<html>Enter the classifier name ' ...
                    '<br>(e.g. EPP, FlowSOM, MLP, etc.)</html>']);
                bp=Gui.BorderPanel([],0,0, 'West', west, 'East', east);
            end
            
            function [jtf, pnl]=make(name, value, tip)
                jtf=Gui.NewTextField(value, 8, tip, ...
                    [],[], @(txt, jtf, e)notEmpty(txt, name), [], true);
                if isempty(char(jtf.getText.trim))
                    jtf.setForeground(Gui.ERROR_COLOR);
                    jtf.setBackground(Gui.WARNING_COLOR);
                end
                pnl=Gui.Panel(jtf);
                Gui.SetTitledBorder([name ' name ...'], pnl);
                pnl.setToolTipText(tip);
            end

            function complaint=notEmpty(txt, name)
                if isempty(strtrim(txt))
                    complaint=['Please enter a ' name ' value'];
                else
                    complaint=[];
                end
            end
        end

        function [cancelled, file, dataset, classifier, view]=Ask(dataset, classifier, file)
            cancelled=false;
            bp=[];
            view=false;
            if nargin<3
                file=[];
                if nargin<2
                    classifier='';
                    if nargin<1
                        dataset='OMIP-077'; %our new favourite
                    end
                    [bp, jtfDataset, jtfClassifier]=ClassificationTable.EntryPanel(dataset);
                end
            end
            readOnly=nargin==1;
            if isempty(bp)
                [bp, jtfDataset, jtfClassifier]=ClassificationTable.EntryPanel(dataset, classifier);
            end
            app=BasicMap.Global;
            if isempty(file)
                priorFile=app.get(ClassificationTable.PROP_FOLDER, ...
                    File.Documents);
                if exist(priorFile, 'dir')
                    priorFile=fullfile(priorFile, app.get(...
                        ClassificationTable.PROP_FILE, ...
                        '5_classifiers.xls'));
                end
            else
                priorFile=file;
            end
            btn=Gui.NewBtn(Html.WrapSmall('Colors'), ...
                @(h,e)ClassificationTable.AskToDownloadColors(), ...
                'Download cell type colors for demos', 'colorWheel16.png');
            chbViewToo=Gui.CheckBox(Html.WrapSmallBold('See plots of all matches NOW.'), ...
                true, [], [], [], ['Show plots for ALL matches dones ' ...
                'on this dataset']);
            if ~readOnly
                myLastFile=['<html>Add into my <i>last used</i> file:' ...
                    Html.FileTree(priorFile) '<hr></html>'];
                ttl='Add the current QFMatch into which file ??';
                choices={...
                    '<html><hr>Download our demo file and add to it</html>', ...
                    '<html><hr>Download our 2nd demo file with EPP included and add to it</html>', ...
                    '<html><hr>I will choose a file<hr></html>', ...
                    myLastFile};
                btnPnl=Gui.FlowLeftPanel(0, 0, btn, chbViewToo);
            else
                myLastFile=['<html>My <i>last</i> file:' ...
                    Html.FileTree(priorFile) '<hr></html>'];
                ttl='Use the "match collection" file from';
                choices={...
                    '<html><hr>Google Cloud (our demo file)</html>', ...
                    '<html><hr>Google Cloud (our 2nd demo file with EPP)</html>', ...
                    '<html><hr>My local file system <hr></html>', ...
                    myLastFile};
                btnPnl=Gui.Panel(btn);
            end
            bp.add(btnPnl, 'South');
            if exist(priorFile, 'file')
                dflt=4;
            else
                choices(end)=[];
                dflt=3;
            end
            while true
                MatBasics.RunLater(@(h,e)focus(), .2);
                [answer, cancelled]=Gui.Ask( ...
                    ttl, choices, ...
                    '', 'Confirm...', dflt, bp);
                if cancelled
                    return;
                end
                if answer==0
                    msgWarning('No choice made');
                end
                dataset=char(jtfDataset.getText.trim);
                if ~isempty(jtfClassifier) 
                    classifier= char(jtfClassifier.getText.trim);
                end
                if isempty(dataset)
                    msgWarning('Enter the dataset name', 5, 'south+');
                elseif ~isempty(jtfClassifier) && isempty(classifier)
                    msgWarning('Enter the classifier name', 5, 'south+');
                else
                    break;
                end

            end
            if answer==4
                file=priorFile;
            elseif answer==1 || answer==2
                if answer==1
                    nClassifiers=5;
                else
                    nClassifiers=6;
                end
                file=ClassificationTable.DownloadClasses(nClassifiers, true);
                if isempty(file)
                    cancelled=true;
                    return;
                end
                [p, f ,e]=fileparts(file);
                app.set(ClassificationTable.PROP_FOLDER, p);
                app.set(ClassificationTable.PROP_FILE, [f e]);                
            else
                file=[];
            end
            view=chbViewToo.isSelected;

            function focus
                jtfDataset.requestFocus;
                jtfDataset.selectAll
                drawnow;
            end
        end

        function AskToDownloadColors
            nEditors=length(BasicMap.Global.cbn.listeners);
            if nEditors>0
                str=Html.WrapSmallBoldOnly(['<br>(<font color="red">'...
                    String.Pluralize2('color editor window', nEditors) ...
                    ' </font>will need to be closed/re-opened);']);
            else
                str='';
            end
            colorsToo=...
                askYesOrNo(Html.SprintfHr(['Download latest ' ...
                'cell type colors?<br><br>(%s)'], ...
                Html.WrapSmallBoldOnly(['This will <font color="red">' ...
                'overwrite</font> your current cell type colors!' str])), ...
                'Get colors?', 'center', false, [], ...
                'ClassificationTable.AskToDownloadColors');
            if colorsToo
                ClassificationTable.DownloadColors(true);
            end
        end

    end
end