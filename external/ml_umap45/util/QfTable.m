    classdef QfTable < handle
%   Class for Hi-D matching with merging of data subsets using QFmatch or
%   F-measure or both.
%
%   This produces a visual table and histogram of dissimilarities and
%   overlaps between subset matches.
%
%   The function getAverages answers the question of overall goodness by
%   returning the median/mean dissimilarity or overlap

%   QF Algorithm is described in 
%   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6586874/
%   and
%   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5818510/
%
%   Bioinformatics lead: Darya Orlova <dyorlova@gmail.com>
%   Software Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties
        contextDescription;
        context;
        fcnFpnSubsets;
        fncSelect;
        btns;
        btnLbls;
    end
    
    properties(SetAccess=private)
        fig;
        data;
        unmatched;
        matchesStr;
        sortTable;
        tb;
        qf;
        tClrs;
        priorFig;
        f1Fig;
        ijiFig;
        rcFig
        similarityFig;
        cFig;
        fHistFig;
        qHistFig;
        rcHistFig;
        cHistFig;
        mcHistFig;
        falsePosNegFig;
        confusionFig;
        otherFigs={};
        R;
        externalArgs;
        app;
        histResize;
        predictions;
        predictionsOfThese;
        predictionsTable;
        listener;
        cbSyncKld;
        cbStackedPredictions;
        cbFlashlight;
        nameIdx;
        szIdx;
        freqIdx;
        symIdx;
        similarityIdx;
        overlapIdx;
        rankIdx;
        concordanceIdx;
        mahalanobisConcordanceIdx
        robustConcordanceIdx;
        idIdx;
        groupIdx;
        fncPredictionSelected;
        gt;
        jdHeatMap;
        idxsHeatMap;
        idMap;
        matchedNames;
        matchedClrs;
        matchMap;
        summaryT;
        saveFile;
        overlapJd;
    end
    
    properties(Constant)
        PROP_OUTER_POS2='OuterPosition2_V2';
        DEBUG=false;
        FIG_NAME='QFMatch';
        FIG_NAME_PRED='PredictionAdjudicator';
    end
    
    methods
        function this=QfTable(qf, tClrs, props, priorFig, ...
                visible, externalArgs, propSuffix, predictions)
            if nargin<8
                predictions=[];
                if nargin<7
                    propSuffix='';
                    if nargin<6
                        externalArgs=[];
                        if nargin<5
                            visible=true;
                            if nargin<4
                                priorFig=get(0, 'currentFig');
                                if nargin<3
                                    props=[];
                                end
                            end
                        end
                    end
                elseif ~startsWith(propSuffix, '.')
                    propSuffix=['.' propSuffix];
                end
            end
            try
                qf.args.pu.setText('Setting up the QFMatch table ...');
            catch
            end
            this.predictions=predictions;
            isForPrediction=~isempty(predictions);
            if isempty(propSuffix)
                if isForPrediction
                    propSuffix='.PredictionsV2';
                end
            end
            if iscell(visible)
                locate_fig=visible;
                visible=true;
            else
                locate_fig={};
            end
            
            this.externalArgs=externalArgs;
            app=BasicMap.Global;
            this.app=app;
            if app.highDef
                this.histResize=.76;
            else
                this.histResize=.66;
            end
            if isempty(props)
                props=app;
            end
            this.tClrs=tClrs;
            this.qf=qf;
            PROP_COL_W=['qftColumnWidthsV2' propSuffix];
            PROP_COL_ORD=['qftColumnOrderV1' propSuffix '_' ...
                num2str(qf.fMeasuringUnmerged)];
            PROP_SORT=['qftRowOrderV1' propSuffix];
            path=BasicMap.Path;
            if isForPrediction
                tableName=QfTable.FIG_NAME_PRED;
            else
                tableName=QfTable.FIG_NAME;
            end
            if visible
                if ~isempty(priorFig) && isempty(locate_fig)
                    pu=PopUp(['Preparing ' tableName ' table'],...
                        'north west+',  'Note...', false, true, ...
                        Gui.GetResizedImageFile('orlova.png', ...
                        .25, BasicMap.Global));
                else
                    pu=[];
                end
            else
                pu=[];
            end
            this.priorFig=priorFig;
            [this.fig, tb_]=Gui.Figure;
            this.tb=tb_;
            if strcmpi(qf.args.mergers, 'test set')
                n2n='1-to-many';
            elseif strcmpi(qf.args.mergers, 'training set')
                n2n='many-to-1';
            elseif strcmpi(qf.args.mergers, 'both')
                n2n='many-to-many';
            else
                n2n='1-to-1';
            end
            figName=[tableName ' ' n2n ' '...
                num2str(length(qf.tIds)) ' X ' num2str(length(qf.sIds))...
                ' subsets'];
            if isfield(this.qf.args, 'windowTitleSuffix') && ~isempty(this.qf.args.windowTitleSuffix)
                figName=[figName ' ' this.qf.args.windowTitleSuffix];
            end
            set(this.fig, 'CloseRequestFcn', @(h, e)hush(h), ...
                'Name', figName);
            this.syncProperties;

            [this.data, labels, fmts, tips,  this.unmatched, groupIdx, ...
                freqIdx, rankIdx, symIdx, this.matchesStr, numIdx,...
                nameIdx, matchesIdx, similarityIdx, overlapIdx, ...
                idIdx, szIdx, concIdx, mahalanobisConcIdx, ...
                robustConcIdx, this.idMap]...
                =QfTable.Contents(qf, tClrs, pu);
            this.nameIdx=nameIdx;
            this.szIdx=szIdx;
            this.freqIdx=freqIdx;
            this.symIdx=symIdx;
            this.similarityIdx=similarityIdx;
            this.overlapIdx=overlapIdx;
            this.rankIdx=rankIdx;
            this.groupIdx=groupIdx;
            this.concordanceIdx=concIdx;
            this.mahalanobisConcordanceIdx=mahalanobisConcIdx;
            this.robustConcordanceIdx=robustConcIdx;
            this.idIdx=idIdx;
            if isForPrediction
                labels{groupIdx}='<html>Subset<br>type</html>';
                labels{similarityIdx}=['<html>' this.app.smallStart ...
                    'Similarity to<br><b>True Class</b>' this.app.smallEnd '</html>'];
                fmts(similarityIdx)=6;
                fmts(groupIdx,:)=[9 nan];
                fmts(idIdx,:)=[8 nan];
                fmts(nameIdx,:)=[26 nan];
                this.adjustForPredictions(predictions, rankIdx,  groupIdx,...
                    similarityIdx, overlapIdx, idIdx, nameIdx, szIdx, matchesIdx);
            end
            this.R=size(this.data);
            [sData, widths]=SortTable.ToSortableAlignedHtml(this.data, fmts);
            this.sortTable=SortTable(this.fig, sData, ...
                labels, [], @(h,e)select(e), tips);
            st=this.sortTable;
            jt=st.jtable;
            st.setSelectionBar;
            if ismac
                st.uit.FontSize=13;
            end
            N=length(widths);
            for i=1:N
                if app.highDef
                    factor=app.toolBarFactor;
                else
                    factor=1;
                end
                st.setColumnWidth(i, widths(i)*factor)
                st.setColumnWidth(i, widths(i)*factor)
            end
            preSorted=SortTable.SetRowOrder(jt, ...
                BasicMap.GetNumbers(props, PROP_SORT));
            if ~preSorted
                if ~isForPrediction
                    jt.sortColumn(rankIdx-1, true, true);
                    jt.sortColumn(groupIdx-1, false, true);
                    jt.sortColumn(freqIdx-1, false, false);
                else
                    jt.sortColumn(nameIdx-1, true, true); %names
                    jt.sortColumn(groupIdx-1, false, false);
                    jt.sortColumn(freqIdx-1, false, false);
                end
            end
            ToolBarMethods.addButton(tb_, ...
                fullfile(path, 'world_16.png'), ...
                ['<html>See matches and histograms<br>'...
                '(if showing) in default browser</html>'], ...
                @(h,e)browse(this, qf.matrixHtml));
            ToolBarMethods.addButton(tb_, ...
                fullfile(path, 'save16.gif'), ...
                'Save table state to excel or text file', ...
                @(h,e)saveSortedTable(this));
            startingRowOrder=SortTable.GetRowOrder(jt);
            bb=ToolBarMethods.addButton(tb_,fullfile(path, 'table.gif'), ...
                'Restore default row and column order', ...
                @(h,e)defaultOrder());
            ToolBarMethods.addButton(tb_,fullfile(path, 'leftArrowNarrow.png'), ...
                'Shift sorted columns to the left', ...
                @(h,e)shiftSortLeft());
            pnlCb=Gui.FlowLeftPanel(2,0);
                
            img=Html.ImgXy('pseudoBarHi.png', [], .819);
            if isForPrediction
                pnlCb.add(javax.swing.JLabel('  '));
                img2=Html.ImgXy('pinFlashlightTransparent.png', [], .92);
                prop4='qfTable.Flashlight';
                this.cbFlashlight=Gui.CheckBox(...
                    ['<html>' img2 '</html>'], ...
                    app.is(prop4, true), app, prop4, [], ...
                    ['<html>See selections highlighted '...
                    'in other plots.</html>']);
                pnlCb.add(this.cbFlashlight);
                pnlCb.add(javax.swing.JLabel('  '));
                prop1='qfTable.Sync1';
                prop4='qfTable.Sync4';
                this.cbStackedPredictions=Gui.CheckBox(...
                    Html.WrapSmallBold(['4 ' img], app), ...
                    app.is(prop4, true), ...
                    [], '', @(h,e)syncKld([], prop4), ...
                    ['<html>Keep seeing DimensionStacker ' img ...
                    'stacked with <br>4 prediction subsets for 1st '...
                    'selection in <u>this</u> table.</html>']);
                pnlCb.add(this.cbStackedPredictions);
                syncWord=['1 ' img];
                syncDflt=false;
                jt.removeColumn(jt.getColumn(labels{matchesIdx}))
                jt.removeColumn(jt.getColumn(labels{rankIdx}))
                groupIdx=5; 
                freqIdx=6; 
                symIdx=8;
                similarityIdx=3; 
                overlapIdx=4; 
                szIdx=7;
            else
                prop1='qfTable.SyncKld';
                syncWord=['Sync ' img];
                syncDflt=true;
            end
            syncDflt=app.is(prop1, syncDflt);

            ToolBarMethods.addSpeakerButton(tb_);
            this.cbSyncKld=Gui.CheckBox(...
                Html.WrapSmallBold(syncWord, app), ...
                syncDflt, [], '', @(h,e)syncKld(prop1, []), ...
                ['<html>Keep seeing ParameterExplorer ' img ...
                ' for 1st selection in <u>this</u> table</html>']);
            pnlCb.add(this.cbSyncKld);
            ToolBarMethods.addComponent(tb_, pnlCb);
            ToolBarMethods.addSeparator(tb_);
            savedOrder=BasicMap.GetNumbers(props, PROP_COL_ORD);
            if ~isempty(savedOrder)
                try
                    SortTable.SetColumnOrder(jt, savedOrder,...
                        BasicMap.GetNumbers(props, PROP_COL_W));
                catch
                end
            else
                defaultColumnOrder;
            end
            PROP_OUTER_POS=[QfTable.PROP_OUTER_POS2 propSuffix];
            op=BasicMap.GetNumbers(props, PROP_OUTER_POS);
            if length(op)==4
                newFigPos=Gui.RepositionOnSameScreenIfRequired(op);
                set(this.fig, 'OuterPosition', newFigPos);
                Gui.FitFigToScreen(this.fig);
            elseif ~isempty(priorFig)
                Gui.SetToRight(this.priorFig, this.fig, false, 80);
            end
            if visible
                if ~isempty(pu)
                    pu.close;
                end
                if ~isempty(locate_fig)
                    SuhWindow.Follow(this.fig, locate_fig);
                    SuhWindow.SetFigVisible(this.fig);
                else
                    Gui.SetFigVisible(this.fig);
                end
                drawnow;
                if isForPrediction
                    Gui.ShowBusy(this.fig, ...
                        Gui.YellowSmall('Prediction results<br>true+/false+/false-'),...
                        'wayneMoore2.png', .61, false, 'South', 2);
                else
                    Gui.ShowBusy(this.fig, ...
                        Gui.YellowSmall('Darya''s QFMatch results'),...
                        'orlova.png', .33, false, 'South', 2);
                end
            end
            if preSorted
                jt.unsort;
                SortTable.SetRowOrder(jt, ...
                    BasicMap.GetNumbers(props, PROP_SORT));
            elseif isempty(savedOrder)
                bb.doClick;
            end
            drawnow;
            rowHeight=floor(jt.getRowHeight*1.33);
            set(this.fig, 'ResizeFcn', @(h,e)resizeQft())
            MatBasics.DoLater(@(h,e)heighten(), .31);
            this.qf.tClrs=this.tClrs;
            if ~isForPrediction
                try
                    ToolBarMethods.addButton(tb_, fullfile(path, 'histQF.png'), ...
                        'See QFMatch histogram (mass+distance similarity)', ...
                        @(h,e)doHistQF(this));
                    ToolBarMethods.addButton(tb_, fullfile(path, 'histF.png'), ...
                        'See overlap histogram based on F-measure', ...
                        @(h,e)doHistF(this, true));
                    if isempty(this.qf.falsePosNegs)
                        this.qf.getFalsePosNegRecords;
                    end
                    ToolBarMethods.addSeparator(tb_);
                    ToolBarMethods.addButton(tb_, fullfile(path, ...
                        'phenogram.png'), ...
                        ['<html>HiD subset views for selections:<ul>'...
                        '<li><u>HeatMap</u> ' this.app.supStart ...
                        '(fast earth-mover''s distance).' this.app.supEnd ...
                        '<li><u>Phenograms/QF-tree</u>' this.app.supStart...
                        '(fast earth-mover''s distance).' this.app.supEnd ...
                        '<li><u>MDS</u>' this.app.supStart ...
                        '(multi-dimensional scaling)' this.app.supEnd '</html>'],...
                        @(h,e)seePlots(this, h))
                    if ~isempty(this.qf.falsePosNegs)
                        try
                            this.qf.getFalsePosNegsTab(false, externalArgs);
                        catch ex
                            %this.qf is a structure and not an instance
                            %so the 2 getFalse*() methods were called and
                            %results are stored in structure
                            if ~isstruct(this.qf)
                                ex.getReport
                            end
                        end
                        ToolBarMethods.addButton(tb_, ...
                            fullfile(path, 'plusMinus.png'), ...
                            'See false positive/negative info', ...
                            @(h,e)seePlotsForSameData(this, h));
                        QfHiDM.ReadingsToolBarBtn(tb_, app);
                        html=QfHiDM.ConvertMatchesTex2Html(this.matchesStr);
                        ToolBarMethods.addComponent(tb_, ...
                            javax.swing.JLabel(html))                        
                    end
                catch ex
                    if ~isstruct(this.qf)
                        ex.getReport
                    else
                        this.qf.falsePosNegs=[];
                    end
                end
            else
                b1='<font color="blue">';b2='</font>';
                [testSetWins, nPredicted, means]=this.getPredictionSummary;
                txt=sprintf('False + <i>MORE similar</i> in %s%d/%d%s cases', ...
                    b1, testSetWins, nPredicted, b2);
                tip=sprintf([...
                    'In %s%d/%d%s cases the <b>false +</b> data has ' ...
                    '<i>greater</i> mass + distance similarity to the '...
                    '<b>True Class</b> than does the <b>false -</b>' ...
                    ' data</td></tr><tr><td>'...
                    '<b>false +</b> are assignments '...
                    'made <i>only</i> by the 2nd classifier (e.g. test ' ...
                    'set, often automatic gating).<br>' ...
                    '<b>false -</b> are assignments made <i>only</i> '...
                    'by the 1st classifier (e.g. training set, ' ...
                    'often manual gating).'...
                    '</td></tr><tr><td>Mean mass + distance similarity '...
                    'of <b>true +</b> / <b>false  +</b> / <b>false -'...
                    '</> is  %s%3.1f%% / %3.1f%% / %3.1f%%%s'], ...
                    b1, testSetWins, nPredicted, b2, ...
                    b1, means(1), means(2), means(3), b2 );
                tip=Html.WrapTable(tip, 3, 5.3, '1', 'left', 'in', app);
                
                [~, jl]=Gui.ImageLabel(Html.WrapSmallBold(txt, app), ...
                    'smallGenie.png', tip, @describeMoreSimilar);
                ToolBarMethods.addComponent(tb_, jl);
            end
            if isfield(this.qf.args, 'highlighter_registry') &&  ~isempty(this.qf.args.highlighter_registry)
                feval(this.qf.args.highlighter_registry, ...
                    @(g,o,r)hearHighlighting(g,o, r));
            end
            Gui.BlowWindow(tb_.jToolbar);

            function ok=hearHighlighting(gate, ~, ~)
                if ~ishandle(this.fig)
                    ok=false;
                    return;
                end
                ok=true;
                dataIdx=this.idMap.get([gate.id(8:end) '.idx']);
                if isempty(dataIdx)
                    return;
                end
                sym=this.data{dataIdx, this.symIdx};
                newSym=regexprep(sym, 'color="#(.*)"', Html.HexColor(gate.highlightColor));
                this.sortTable.setValue(dataIdx, this.symIdx, newSym);
                this.idMap.set([gate.id(8:end) '.clr'], gate.highlightColor);
                this.data{dataIdx, this.symIdx}=newSym;
            end

            function describeMoreSimilar(h,~)
                jw=Gui.WindowAncestor(h);
                msg(struct( 'javaWindow', jw,'msg', tip, ...
                    'icon', 'genieSearch.png'), ...
                    14, 'north east+', 'Possibly useful info...');
            end

            function syncKld(prop1, prop4)
                go1=false; go4=false;
                if ~isempty(prop1)
                    if this.cbSyncKld.isSelected
                        go1=true;
                        app.set(prop1, 'true');
                    else
                        app.set(prop1, 'false');
                    end
                end
                if ~isempty(prop4)
                    if this.cbStackedPredictions.isSelected
                        go4=true;
                        app.set(prop4, 'true');
                    else
                        app.set(prop4, 'false');
                    end
                end
                if go1 || go4                    
                    try
                        this.listener.reselect;
                    catch ex
                        ex.message
                    end
                end
            end

            function resizeQft()
                MatBasics.DoLater(@(h,e)budge(), .31);
                function budge
                    jt.setRowHeight(rowHeight);
                    drawnow;
                end
            end

            function heighten
                drawnow;
                jt.setRowHeight(rowHeight);
                try
                    if props.multiProps.is('mds_figHistF', false)
                        this.doHistF;
                    end
                    if props.multiProps.is('mds_figHistQF', false)
                        this.doHistQF(visible);
                    end
                catch ex
                    %disp('Using QF table outside of CytoGenie''s AutoGate')
                    %ex.getReport
                end
            end
            
            function shiftSortLeft()
                if ~SortTable.MoveSortColumnsLeft(jt)
                    if ismac
                        key='command';
                    else
                        key='control';
                    end
                    msg(Html.Wrap(['There is no sort order currently.'...
                        '<br>To sort the table you<ul>'...
                        '<li>Click on a column to sort it'...
                        '<li>Click again for descending order'...
                        '<li>Hold ' key ' key for multi column</ul>']), ...
                        8, 'north', 'No sort order...');
                else
                    r=jt.getSelectedRow;
                    if r<0
                        r=0;
                    end
                    rect=jt.getCellRect(r,0, true);
                    jt.scrollRectToVisible(rect);
                end
            end
            
            function defaultOrder()
                defaultRowOrder;
                defaultColumnOrder;
            end
            
            function defaultRowOrder
                if isForPrediction
                   %NOT -1 for groupidx and freqIdx since we removed 
                   % matches column for predictions
                    SortTable.SetRowOrder(jt, [(nameIdx-1) true 0; ...
                        groupIdx false 0; freqIdx false 1]);
                else
                    SortTable.SetRowOrder(jt, [(rankIdx-1) true 0; ...
                        (groupIdx-1) false 0; (freqIdx-1) false 1]);
                end
            end
            
            function defaultColumnOrder()
                idxs=1:jt.getColumnCount;
                if isForPrediction
                    eliminate=[nameIdx, overlapIdx, similarityIdx,...
                        numIdx,  symIdx, szIdx, freqIdx];
                    %+1 since we removed matches column for predictions
                    startWith=[nameIdx, similarityIdx+1, szIdx+1, ...
                        freqIdx+1, overlapIdx+1];
                else
                    if ~this.qf.fMeasuringUnmerged
                        startWith=[rankIdx symIdx numIdx nameIdx similarityIdx matchesIdx groupIdx freqIdx szIdx];
                    else
                        startWith=[rankIdx symIdx numIdx nameIdx overlapIdx];
                    end
                    eliminate=startWith;                        
                end
                eliminate=sort(eliminate, 'descend');
                idxs(eliminate)=[];
                idxs=[startWith idxs]-1;
                SortTable.SetColumnOrder(jt, idxs, [], true);
            end
            
            
            function hush(h)
                try
                    [columnOrder, widths]=SortTable.GetColumnOrder(jt);
                    props.set(PROP_COL_ORD, num2str(columnOrder));
                    props.set(PROP_COL_W, num2str(widths));
                    rowOrder=SortTable.GetRowOrder(jt);
                    if ~isequal(rowOrder, startingRowOrder)
                        props.set(PROP_SORT, MatBasics.Encode(...
                            rowOrder));
                    end
                    if exist('PROP_OUTER_POS', 'var') %r2017a thread issue
                        props.set(PROP_OUTER_POS, ...
                            num2str(get(h, 'OuterPosition')));
                    end
                catch 
                end
                if ~isempty(this.priorFig) && ishandle(this.priorFig)
                    figure(this.priorFig);
                    Gui.CloseFigs(this.otherFigs);
                end
                delete(h);
            end
            
            function select(e)
                if e.getValueIsAdjusting || isempty(this.fncSelect)
                    return
                end
                colIdx=st.jtable.getSelectedColumn;
                rowIdxs=st.jtable.getSelectedRows;
                if ~isempty(rowIdxs)
                    try
                        colIdx=jt.convertColumnIndexToModel(colIdx);
                        if isempty(this.predictions)
                            st.showTip(colIdx+1);
                        end
                        N_=length(rowIdxs);
                        isTeachers=zeros(1,N_);
                        qfIdxs=zeros(1,N_);
                        nTid=length(this.qf.tIds);
                        for j=1:N_
                            rowIdx=rowIdxs(j);
                            rowIdx=jt.getActualRowAt(...
                                jt.convertRowIndexToModel(rowIdx))+1;
                            isTeachers(j)=rowIdx<=nTid;
                            if isTeachers(j)
                                qfIdxs(j)=rowIdx;
                            else
                                qfIdxs(j)=rowIdx-nTid;
                            end
                        end
                        try
                            try
                                feval(this.fncSelect, this.qf, ...
                                    isTeachers, qfIdxs);
                            catch ex
                                %backward compatability
                                feval(this.fncSelect, this.qf, ...
                                    isTeachers(1), qfIdxs(1));
                            end
                        catch ex
                            ex.getReport
                            if isTeachers(1)
                                name=this.qf.tNames{qfIdxs(1)};
                            else
                                name=this.qf.sNames{qfIdxs(1)};
                            end
                            fprintf('rowIdx=%d, teacher=%d, name="%s"\n',...);
                                rowIdx, isTeachers(1), name);
                        end
                    catch ex
                        ex.getReport
                    end
                    return;
                else
                    try
                        feval(this.fncSelect, this.qf, ...
                            [], []);
                    catch ex
                        ex.getReport
                    end

                end
                st.showTip(0);
            end
        end
        
        function syncProperties(this)
            try
                args=this.qf.args;
                isClass=~isstruct(this.qf);
                if isa(args.properties, 'MultiProps')
                    % find non GLOBAL properties
                    props=args.properties.instances{2};
                else
                    props=args.properties;
                end
                go('dataSetName', args.dataSetProperty);
                go('classifierName', args.classifierProperty);
                if ~isempty(this.qf.dataSetName)
                    name=get(this.fig, 'Name');
                    if contains(name, ', dataset: ')
                        idx=String.IndexOf(name, ', dataset: ');
                        set(this.fig, 'Name', [name(1:idx-1) ', dataset: ' ...
                            this.qf.dataSetName ', classifier: ' ...
                            this.qf.classifierName]);
                    else
                        set(this.fig, 'Name', [name ', dataset: ' ...
                            this.qf.dataSetName ', classifier: ' ...
                            this.qf.classifierName]);
                    end
                end
            catch ex
                ex.getReport
            end

            function go(field, property)
                if ~isempty(property)
                    if ~isClass
                        fieldValue=getfield(this.qf, field);
                    else
                        fieldValue=this.qf.(field);
                    end
                    if isempty(fieldValue)
                        propertyValue=props.get(property);
                        if ~isempty(propertyValue)
                            if isClass
                                this.qf.(field)=propertyValue;
                            else
                                this.qf=setfield(this.qf, field, propertyValue);
                            end
                        end
                    else
                        props.set(property, fieldValue);
                    end
                end
            end
        end

        function t=getTableData(this)
            t=this.sortTable.getTableData(this.data);
        end

        function setPredictionListener(this,fnc)
            this.fncPredictionSelected=fnc;
        end
        
        function adjustForPredictions(this, predictions, rankIdx, plotIdx,...
                    similarityIdx, overlapIdx,  idIdx, nameIdx, szIdx, matchesIdx)
            R_=size(this.data);
            nT=predictions.nTeachers;
            sfxPre=[this.app.supStart '&nbsp;&nbsp;<b><i>True Class</i></b>' this.app.supEnd '</html>'];
            sfxFn=['</font>' this.app.supStart ' <b>false <font color="red">-</font></b>' this.app.supEnd '</html>'];
            sfxFp=['</font>' this.app.supStart ' <b>false <font color="red">+</b></font>' this.app.supEnd '</html>'];
            sfxTp=['</font>' this.app.supStart ' <b>true +</b>' this.app.supEnd '</html>'];
            for i=1:nT
                this.data{i, plotIdx}='predicted';
                this.data{i, rankIdx}=nan;
                this.data{i, similarityIdx}=nan;
                this.data{i, overlapIdx}=nan;
                this.data{i, matchesIdx}=nan;
                name=this.data{i,nameIdx};
                this.data{i, nameIdx}=['<html><' ...
                    strrep(strrep(strrep(name, '<',''),'>',''), '/','') ...
                    ' >' name sfxPre];
            end
            nPredictions=R_-nT;
            if SuhPredictions.DEBUG
                assert(nPredictions==length(predictions.sNames));
            end
            predictions.compress;
            for i=1:nPredictions
                i2=nT+i;
                this.data{i2, rankIdx}=nan;
                name=this.data{i2, nameIdx};
                strId=this.data{i2, idIdx};
                id=str2double(strId);
                [similarity, overlap]=predictions.describe(id);
                this.data{i2, similarityIdx}=similarity;
                this.data{i2, overlapIdx}=overlap;
                if endsWith(strId, '.3')
                    word='false -';
                    name=['<html><' ...
                        strrep(strrep(strrep(name, '<',''),'>',''), '/','') ...
                        ' >&nbsp;&nbsp;<font color="#B66666">' name(1:end-7) sfxFn];
                elseif endsWith(strId, '.2')
                    word='false +';
                    name=['<html><' ...
                        strrep(strrep(strrep(name, '<',''),'>',''), '/','') ...
                        ' >&nbsp;&nbsp;<font color="#B66666">' name(1:end-7) sfxFp];
                else
                    word='true +';
                    name=['<html><' ...
                        strrep(strrep(strrep(name, '<',''),'>',''), '/','') ...
                        ' >&nbsp;&nbsp;<font color="#66B666">' name(1:end-6) sfxTp];
                end
                this.data{i2, nameIdx}=name;
                this.data{i2, plotIdx}=word;
                this.data{i2, matchesIdx}=nan;
                if SuhPredictions.DEBUG
                    if endsWith(strId, '.3')
                        sizeExact=sum(predictions.negLbls==id);
                    else
                        sizeExact=sum(predictions.posLbls==id);
                    end
                    c=int32(2+(id-floor(id))*10);
                    r=predictions.sums(:,1)==floor(id);
                    assert(sizeExact==predictions.sums(r,c))
                    sizeRough=this.data{i2,szIdx};
                    similarityRough=this.data{i2,similarityIdx};
                    overlapRough=this.data{i2,overlapIdx};
                    if sizeExact ~= sizeRough
                        fprintf('"%s" size, from %d to %d\n', ...
                            predictions.sNames{i}, sizeRough, sizeExact)
                    end
                    if similarityRough ~= similarity
                        fprintf('"%s" size, from %d to %d\n', ...
                            predictions.sNames{i}, similarityRough, ...
                            similarity)
                    end
                    if overlapRough ~= overlap
                        fprintf('"%s" size, from %d to %d\n', ...
                            predictions.sNames{i}, overlapRough, ...
                            overlap)
                    end
                end
            end
            predictions.decompress;
        end
            
        function listener=listen(this, columnNames, ...
                trainingSet, testSet, trainingIds, testIds, ...
                trainingSetName, testSetName, parseFileNameOut)
            if nargin>4
                if nargin<9
                    parseFileNameOut=false;
                    if nargin<8
                        testSetName='';
                        if nargin<7
                            trainingSetName='';
                        end
                    end
                end
                listener=SuhMatchListener(this, columnNames, ...
                    trainingSet, testSet, trainingIds, testIds, ...
                    trainingSetName, testSetName, parseFileNameOut);
            else
                if nargin<4
                    listener=SuhMatchListener(this, ...
                        columnNames, trainingSet);
                else
                    listener=SuhMatchListener(this, columnNames, ...
                        trainingSet, testSet);
                end
            end
            this.fncSelect=@(qf,isTeachers, idxs)select(...
                listener, qf,isTeachers,idxs);
            if ~isempty(this.qf) 
                if ~isstruct(this.qf) || isfield(this.qf, 'columnNames')
                    if isempty(this.qf.columnNames)
                        this.qf.columnNames=columnNames;
                    end
                end
            end
            this.listener=listener;
        end
        
        function [medianSimilarity, meanSimilarity, ...
                medianOverlap, meanOverlap, ...
                medianConcordance, meanConcordance,...
                medianRobustConcordance, meanRobustConcordance]...
                =getAverages(this)
            [~,medianSimilarity, meanSimilarity]=this.getData(1);
            [~,medianOverlap, meanOverlap]=this.getData(0);
            [~,medianConcordance, meanConcordance]=this.getData(2);
            [~,medianRobustConcordance, meanRobustConcordance]...
                =this.getData(3);
        end

        function matrix=getMatrix(this, columnIndex)
            if nargin<2
                columnIndex=this.similarityIdx;
            end
            try
                matrix=cell2mat(this.data(:, columnIndex));
            catch ex
                ex.getReport
                matrix=[];
            end
        end
        
        function [testSetWins, nPredicted, means, stddevs, mdns, results]...
                =getPredictionSummary(this)
            if this.isPredictions                
                similarities=this.getMatrix;
                similarities=similarities*100;
                N=length(similarities);
                results=nan(N,4);
                nPredicted=0;
                notPredictedIdxs=[]; notPredictedIdxs(N) = 0;
                nNotPredicted = 0;
                for i=1:N
                    strId=this.data(i,this.idIdx);
                    id=str2double(strId);
                    predictedId=floor(id);
                    idx=find(results(:,1)==predictedId,1);
                    if isempty(idx)
                        nPredicted=nPredicted+1;
                        idx=nPredicted;
                        results(idx,1)=predictedId;
                    end
                    if endsWith(strId, '.3')
                        results(idx, 4)=similarities(i);
                        if similarities(i)==100
                            nNotPredicted = nNotPredicted + 1;
                            notPredictedIdxs(nNotPredicted)=idx;
                        end
                    elseif endsWith(strId, '.2')
                        results(idx, 3)=similarities(i);
                    elseif endsWith(strId, '.1')
                        results(idx, 2)=similarities(i);
                    end
                end
                notPredictedIdxs = notPredictedIdxs(1:nNotPredicted);
                results=results(1:nPredicted,:);
                nums=results(:, 2:4);
                testSetWins=0;
                for i=1:nPredicted
                    %nums(i,:)
                    if isnan(nums(i,1)) && isnan(nums(i,2))
                        % no + prediction true or false then training set wins
                    elseif isnan(nums(i,2)) % if NO false+ but true+ 
                        if nums(i,1)>nums(i,3)%IF true+ more similar than false-?
                            testSetWins=testSetWins+1;
                        end
                    elseif isnan(nums(i,3)) % NO false- then true+ MUST equal training set
                        testSetWins=testSetWins+1;
                    elseif nums(i,2)>=nums(i,3) 
                        % false+ more similar to training set than false-
                        % or tied then the test set wins
                         testSetWins=testSetWins+1;
                    else %false- more similar than false+ training set wins
                    end
                end
                if ~isempty(notPredictedIdxs)
                    nums(notPredictedIdxs,:)=[];
                    %either the training subset had zero matches 
                    % or had another stronger training subset in the same match
                    % group that stole all the matches ...
                    % Hence on Nov 7,2021 I figured out a sharing agreement
                    nPredicted=nPredicted-length(notPredictedIdxs);
                end
            else
                testSetWins=[];
                nPredicted=0;
                nums=[];
                results=[];
            end
            if nargout>2
                nums(isnan(nums))=0;
                means=mean(nums);
                if nargout>3
                    stddevs=std(nums);
                    if nargout>4
                        mdns=median(nums);
                    end
                end
            end
        end        
        
        function [matrix, mdn, mn, sizes]=getData(this, metric)
            if nargin<2 
                metric=1;
            end
            if metric==1
                columnIndex=this.similarityIdx;
            elseif metric==3
                columnIndex=this.mahalanobisConcordanceIdx;
            elseif metric==4
                columnIndex=this.robustConcordanceIdx;
            elseif metric==2
                columnIndex=this.concordanceIdx;
            else
                columnIndex=this.overlapIdx;
            end
            if isnan(columnIndex) || columnIndex<1
                matrix=[];
                mdn=nan;
                mn=nan;
                sizes=nan;
                return;
            end
            matrix=cell2mat(this.data(:, columnIndex))*100;
            matrix(isnan(matrix))=0;
            if any(matrix<0)
                matrix=matrix(matrix>=0);
            end
            trainers=matrix(strcmp('yes',this.data(:, this.groupIdx)));
            % IF multiple trainers are merged to match a test set 
            %    then this will NOT agree with MatBasics.ConfusionF1
            mdn=median(trainers);
            mn=mean(trainers);
            sizes=cell2mat(this.data(:, this.szIdx));
        end
        
        function yes=isPredictions(this)
            yes=~isempty(this.predictions);
        end
        
        function predictionsOfThese=getPredictionsOfThese(this)
            if ~isempty(this.predictionsOfThese)
                predictionsOfThese=this.predictionsOfThese;
            else
                if isstruct(this.qf)
                    if isfield(this.qf, 'predictions')
                        predictionsOfThese=this.qf.predictions;
                        if ~isempty(predictionsOfThese)
                            predictionsOfThese.setMatch(this);
                        end
                    else
                        msg('Re-build the match to see predictions');
                        predictionsOfThese=[];
                    end
                else
                    predictionsOfThese=SuhPredictions(this);
                end
                this.predictionsOfThese=predictionsOfThese;
            end
            if ~isempty(this.fncPredictionSelected)
                this.predictionsOfThese.setSelectionListener(...
                    this.fncPredictionSelected);
            end
        end
        
        function [table, predictions]=seePredictionOfThese(this)
            predictions=this.getPredictionsOfThese;
            if ~isempty(predictions)
                try
                    table=predictions.showTable;
                catch ex
                    table=[];
                    predictions=[];
                    ex.getReport
                    msg(Html.WrapHr(['Re-run the matching ... '...
                        '<br>the current cache is <br>'...
                        'based on an old version!']));
                end
            else
                table=[];
            end
            this.predictionsTable=table;
        end
        
        function showHeatMap(this)
            if isempty(this.qf.tMdns)
                msgError(['<html>Re-run QFMatch with this AutoGate'...
                    '<br>version to get HeatMap data!</html>'], ...
                    8, 'south', 'Prior version made match...');
                return;
            end
            pu=PopUp('Building HeatMap', 'north west');
            jt=this.sortTable.jtable;
            [dataIdxs, names, freqs, syms, mdns, compMdns]=...
                QfHiDM.GetDataInTableOrder(jt, this.qf, this.data, ...
                this.nameIdx, this.freqIdx, this.symIdx, this.idIdx, ...
                this.rankIdx, true); 
            if ~isempty(this.jdHeatMap)
                if ~isequal(this.idxsHeatMap, dataIdxs)
                    this.jdHeatMap.dispose;
                    this.jdHeatMap=[];
                else
                    if ~this.jdHeatMap.isVisible
                        this.jdHeatMap.setVisible(true);
                    end
                    this.jdHeatMap.requestFocus;
                    pu.close;
                    return;
                end
            end
            mNames=this.qf.columnNames;
            N=size(this.qf.tMdns, 2);%avoid classification label 
            for i=1:N
                name=mNames{i};
                idx=String.IndexOf(name, ':');
                if idx>0
                    mNames{i}=name(1:idx-1);
                end
            end
            this.idxsHeatMap=dataIdxs;
            this.jdHeatMap=SuhHeatMap.New(...
                'measurements', mdns, 'rawMeasurements', compMdns,...
                'measurementNames', mNames, ...
                'subsetName', ['Subset ' this.app.supStart 'rank'...
                this.app.supEnd], 'subsetSymbol', syms, ...
                'names', names, 'freqs', freqs, ...
                'cellClickedCallback', @rowClicked,...
                'rowClickedCallback', @rowClicked,...
                'parentFig', this.fig, ...
                'rowClickAdvice', '(<i>click to sync match table</i>)',...
                'windowTitleSuffix', ' for matched subsets');
            SuhWindow.Follow(this.jdHeatMap, this.fig, 'west++');
            pu.close;
            
            function rowClicked(~, row)
                vc=jt.convertColumnIndexToView(idIdx-1);
                dataIdx=dataIdxs(row);
                id=this.data{dataIdx, idIdx};
                R_=jt.getRowCount;
                for r=0:R_-1
                    if strcmp(id, char(jt.getValueAt(r, vc)))
                        jt.scrollRowToVisible(r)
                        jt.setRowSelectionInterval(r,r);
                        figure(this.fig)
                    end
                end
            end
            
        end

        function [name2, s]=get2ndBestMatchName(this, sId, name1, s)
            if isempty(s)
                mf=this.qf.matrixFinal;
                s.mf=mf(1:length(this.qf.sNames), 1:length(this.qf.tNames));
                [s.keyRows, s.keyCols]=MatBasics.Min(s.mf);
            end
            sIdx=find(this.qf.sIds==str2double(sId), 1);
            all=find(s.keyCols==sIdx);
            if length(all)>1
                val=s.mf(sIdx, all);
                [~, I]=sort(val);
                sIdx2nd=all(I(2));
                name2=this.qf.tNames{sIdx2nd};
                if strcmp(name1, name2)
                    sIdx2nd=all(I(1));
                    name2=this.qf.tNames{sIdx2nd};
                    if strcmp(name1, name2)
                        name2='';
                    end
                end
            else
                name2='';
            end
        end
        
        function [names, clrs, map]=getMatchedNames(this, mustOverlap)
            if nargin<2
                mustOverlap=true;
            end
            if isempty(this.matchedNames)
                [names, clrs, sIds]=LabelBasics.NamesAndColors(...
                    this.qf.sIdPerRow(:,1), this.idMap);
                map=Map;
                [matchNames, ~, hasMatch]=...
                    QfHiDM.GetMatchingNamesAndColors(this.qf, ...
                    this.idMap, mustOverlap);
                NN=length(matchNames);
                if NN==length(names)
                    if sIds(1)==0
                        sIdx=1;
                    else
                        sIdx=0;
                    end
                    for i=1:NN
                        mn=matchNames{i};
                        if ~isempty(mn)
                            name=names{i};
                            if startsWith(name, [mn ' (']) ...
                                    && endsWith(name, ')')
                                idx2=length(mn)+3;
                                name(end)='}';
                                names{i}=['\bf' mn ' ^{' name(idx2:end)];
                                %disp('%already done')
                            else
                                names{i}=['\bf' mn ' ^{' names{i} '}'];
                            end
                            if hasMatch(i)
                                map.set(num2str(sIds(sIdx+i)), mn)
                            end
                        end
                    end
                end
                this.matchedNames=names;
                this.matchedClrs=clrs;
                this.matchMap=map;
                if ~isempty(this.qf.falsePosNegs)
                    QfHiDM.GetUnmatchedInOverlapMatrix(map, ...
                        this.qf.matrixOverlap, this.qf.sSizes, this.qf.sIds);
                end
            else
                names=this.matchedNames;
                clrs=this.matchedClrs;
                map=this.matchMap;
            end
        end

        function setPredictionsTable(this,table)
            this.predictionsTable=table;
        end

        function seePlotsForSameData(this, h)
            jMenu=PopUp.Menu;
            assert(this.qf.areEqual);
            scoresMenu=Gui.NewMenu(jMenu, 'Scores X log10(sizes)', ...
                [], 'smallGenie.png');
            Gui.NewMenuItem(scoresMenu, 'Similarity X log10(sizes)', ...
                @(h,e)similarityScoresBySizes(this));
            Gui.NewMenuItem(scoresMenu, 'F1-scores X log10(sizes)', ...
                @(h,e)f1ScoresBySizes(this));
            Gui.NewMenuItem(scoresMenu, 'Concordance X log10(sizes)', ...
                @(h,e)cScoresBySizes(this));
            Gui.NewMenuItem(scoresMenu, ['Inlier concordances X' ...
                ' log10(sizes)'], ...
                @(h,e)ijiScoresBySizes(this));
            Gui.NewMenuItem(scoresMenu, ['Robust inlier concordances ' ...
                'X log10(sizes)'], ...
                @(h,e)rcScoresBySizes(this));
            Gui.NewMenuItem(jMenu, 'False positive/negative charts', ...
                @(h,e)browseFalsePosNeg(this), 'predictions.png');
            Gui.NewMenuItem(jMenu, 'Prediction adjudicator table', ...
                @(h,e)seePredictionOfThese(this), 'table.gif');
            Gui.NewMenuItem(jMenu, 'Confusion chart', ...
                @(h,e)confusionChart(this), 'demoIcon.gif');
            oMnu=Gui.NewMenu(jMenu, 'Overlap matrix for ...', [], 'table.gif');
            Gui.NewMenuItem(oMnu, 'Browser viewing', ...
                @(h,e)showDecisionMatrix(this, true, true));
            Gui.NewMenuItem(oMnu, 'Gating tree interactions', ...
                @(h,e)showDecisionMatrix(this, true, false));
            
            concordanceMenu=Gui.NewMenu(jMenu, ...
                'Concordance (Jaccard index) histograms', ...
                [], 'histF.png');
            Gui.NewMenuItem(concordanceMenu, 'Concordance histogram', ...
                @(h,e)doHistConcordance(this, true));
            Gui.NewMenuItem(concordanceMenu, ['Inlier concordance' ...
                ' histogram'], ...
                @(h,e)doHistMahalanobisConcordance(this, true));
            Gui.NewMenuItem(concordanceMenu, ['Robust inlier ' ...
                'concordance histogram'], ...
                @(h,e)doHistRobustConcordance(this, true));
            Gui.NewMenuItem(jMenu, 'Store this match with other matches...', ...
                @(h,e)integrateClassification(this, true), 'pickGate.png');
                
            jMenu.show(h, 15, 25);
        end

        function seePlots(this, h)
            jMenu=PopUp.Menu;
            jm=Gui.NewMenu(jMenu, ['<html>Phenogram ' this.app.supStart...
                '(QF-tree</i>)' this.app.supEnd '</html>'],...
                [], 'phenogram.png');
            Gui.NewMenuItem(jm, 'Training and test sets', ...
                @(h,e)phenogram(this));
            Gui.NewMenuItem(jm, 'Test set only', ...
                @(h,e)phenogram(this, false));
            Gui.NewMenuItem(jm, 'Training set only', ...
                @(h,e)phenogram(this, true));            
            jm=Gui.NewMenu(jMenu, ['<html>MDS ' this.app.supStart...
                '(multi-dimensional scaling)' this.app.supEnd ...names
                '</html>'], [], 'mds.png');
            hasTestSet= ~isstruct(this.qf) || isfield(this.qf, 'sIdPerRow');
            if ~hasTestSet
                Gui.NewMenuLabel(jm, 'Re-run match for test set')
            end
            mi=Gui.NewMenuItem(jm, 'Training and test sets', ...
                @(h,e)mds(this));
            mi.setEnabled(hasTestSet);
            mi=Gui.NewMenuItem(jm, 'Test set only', ...
                @(h,e)mds(this, false));
            mi.setEnabled(hasTestSet);
            Gui.NewMenuItem(jm, 'Training set', ...
                @(h,e)mds(this, true));
            Gui.NewMenuItem(jMenu, 'Show heat map...', ...
                @(h,e)showHeatMap(this), 'heatMapHot.png');
            oMnu=Gui.NewMenu(jMenu, 'Match score matrix for...', [], 'table.gif');
            Gui.NewMenuItem(oMnu, 'Browser viewing', ...
                @(h,e)showDecisionMatrix(this, false, true));
            Gui.NewMenuItem(oMnu, 'Gating tree interactions', ...
                @(h,e)showDecisionMatrix(this, false, false));
            jMenu.show(h, 15, 25);
        end

        
        function showDecisionMatrix(this, doOverlap, forBrowser, resetOrder)
            if nargin<4
                resetOrder=false;
                if nargin<3
                    forBrowser=false;
                end
            end
            Gui.ShowBusy(this.fig, Gui.YellowH2(...
                'Computing html matrix'), ...
                'Genie.png', 1.1, true);
            if isempty(this.qf.classifierName)
                testTitle='Test set';
                trainingTitle='Training set';
                item='subsets';
            else
                trainingTitle='Reference<br>populations';
                testTitle=[this.qf.classifierName '<br> populations'];
                item='populations';
            end
            if isstruct(this.qf)
                if ~isfield(this.qf, 'overlap_tNames')
                    this.qf.overlap_tNames={};
                    this.qf.overlap_sNames={};
                    this.qf.rowsI=[];
                    this.qf.colsI=[];
                    this.qf.matrixHtmlOverlap=[];
                    this.qf.matrixHtmlOverlap_browser=[];
                end
            end
            try
                [this.qf, html, ttl, generated]=QfHiDM.DecisionMatrix(this.qf, ...
                    false, doOverlap, forBrowser, this.qf.dataSetName,...
                    trainingTitle, testTitle);
                if generated
                    this.save( this.qf, this.saveFile );
                end
                Gui.HideBusy(this.fig);
                if ~isempty(html)
                    if forBrowser
                        Html.BrowseString(Html.Wrap(html));
                    else
                        args=struct(...
                            'fnc', this.fncSelect, 'where', 'west++', ...
                            'title', ttl);
                        args.tNames=this.qf.tNames;
                        args.sNames=this.qf.sNames;
                        args.trainingTtl=strrep(strrep(trainingTitle, ...
                            '<br>', ' '), 'populations', 'population');
                        args.testTtl=strrep(strrep(testTitle, ...
                            '<br>', ''), 'populations', 'population');
                        jd=QfHiDM.TextPane(this.qf, html, args);
                        SuhWindow.Follow(jd, this.fig, 'west++');
                    end
                end
                if doOverlap && (resetOrder || ~Gui.IsVisible(this.overlapJd))
                    try
                        if ~isempty(this.overlapJd)
                            this.overlapJd.dispose;
                        end
                        [this.overlapJd, or]=Gui.Orderer(...
                            this.qf.overlap_tNames, ...
                            this.qf.overlap_sNames, ...
                            ['<html>' strrep(trainingTitle, '<br>', ' ') ...
                            this.app.supStart [' (<b>columns</b>)' ...
                            ''] this.app.supEnd  '</html>'], ...
                            ['<html>' strrep(testTitle, '<br> ',' ') ...
                            this.app.supStart ' (<b>rows</b>)' ...
                            this.app.supEnd  '</html>'],...
                            @(left, right, or)apply(left, right), ...
                            @(oLeft, oRight, or)reset(), ...
                            ['Ordering of overlap matrix ' item], this.fig, ...
                            'south', this.app.adjustHighDef(725, 100), ...
                            this.app.adjustHighDef(300, 170), 'population');
                        or.setIsLeftRightOk(false);
                        MatBasics.RunLater(@(h,e)focus, 4)
                    catch ex
                        ex.getReport
                    end
                end
            catch ex
                Gui.MsgException(ex);
            end
            Gui.HideBusy(this.fig);

            function focus
                this.overlapJd.requestFocus;
            end

            function apply(left, right)
                try
                    this.qf=QfHiDM.ResetOverlapMatrix(this.qf, left, right);
                    this.showDecisionMatrix(doOverlap, forBrowser);
                    this.save( this.qf, this.saveFile, 'overlap order' );
                catch ex
                    Gui.MsgException(ex);
                end
            end

            function [newLeft, newRight]=reset()
                try
                    this.qf=QfHiDM.ResetOverlapMatrix(this.qf);
                    this.showDecisionMatrix(doOverlap, forBrowser, true);
                    this.save( this.qf, this.saveFile, 'overlap order');
                    newLeft=this.qf.overlap_tNames;
                    newRight=this.qf.overlap_sNames;
                catch ex
                    Gui.MsgException(ex);
                end
            end
        end

        function [fig, qf, qft, fig2, qf2, qft2]...
            =phenogram(this, trainingSet, locate)
            if nargin<3
                locate=[];
            end
            if nargin<2
                [fig, qf, qft]=this.phenogram(true);
                [fig2, qf2, qft2]=this.phenogram(false, ...
                    {fig, 'east++', false});
                return;
            end
            if trainingSet 
                [names, clrs]=LabelBasics.NamesAndColors(...
                    this.qf.tIdPerRow(:,1), this.idMap);
                data_=this.qf.tData;
                lbls=this.qf.tIdPerRow(:,1);
                where='south west+';
                ttl='the manually gated training set';                
            else
                [names, clrs]=this.getMatchedNames;
                data_=this.qf.sData;
                lbls=this.qf.sIdPerRow(:,1);
                where='south east+';
                ttl=['the ' this.qf.classifierName ' test set ' ];
            end
            if ~isempty(this.qf.dataSetName)
                ttl={ttl, ['For the ' this.qf.dataSetName ' dataset']};
            else
                ttl={ttl, 'For the dataset *'};
            end
            if isempty(locate)
                locate={this.fig, where, false};
            end
            [qf, qft, fig]=run_QfTree(data_, lbls, ttl,...
                'trainingNames', names, 'log10', true, 'colors', clrs, ...
                'locate_fig', locate);
            %Gui.Add
        end

        function [fig, mds, fig2, mds2]=mds(this, trainingSet, locate)
            if nargin<3
                locate=[];
            end
            if nargin<2
                [fig, mds]=this.mds(true);
                [fig2, mds2]=this.mds(false, ...
                    {fig, 'east++', false});
                return;
            end
            if trainingSet
                [names, clrs]=LabelBasics.NamesAndColors(...
                    this.qf.tIdPerRow(:,1), this.idMap);
                mdns=this.qf.tMdns;
                sizes=this.qf.tSizes;
                where='north west+';
                ttl='MDS: training set populations';
            else
                [names, clrs]=this.getMatchedNames;
                mdns=this.qf.sMdns;
                sizes=this.qf.sSizes;
                where='north east+';
                ttl='MDS: test set populations';
            end
            if isempty(locate)
                locate={this.fig, where, true};
            end
            [mds, fig]=MDS.New(names, this.qf.columnNames, ...
                clrs, mdns, sizes, ttl, locate);
        end

        function [theTitle, hasNames]=getClassificationTitle(this)
            hasNames=true;
            if isempty(this.qf.classifierName)
                figure(this.fig);
                [~, this.qf.classifierName, this.qf.dataSetName]=...
                    ClassificationTable.GetClassifierAndDataSet(...
                    this.qf.classifierName, this.qf.dataSetName, ...
                    'Annotate this match?', 'north');
                this.syncProperties;
            end
            if ~isempty(this.qf.classifierName) &&...
                ~isempty(this.qf.dataSetName)
                theTitle=sprintf('%s run on %s', ...
                    this.qf.classifierName, ...
                    this.qf.dataSetName);
            else
                theTitle='Averages are for trainers (only)';
                hasNames=false;
            end
        end

        function done=f1ScoresBySizes(this, startUp)
            done=true;
            if nargin>1 && startUp
                done=this.scoresBySizes('F1-Score', 'f1Fig', ...
                    this.overlapIdx, 'south++');
            else
                this.scoresBySizes('F1-Score', 'f1Fig', ...
                    this.overlapIdx, 'east+');
            end                        
        end

        function ijiScoresBySizes(this)
            this.scoresBySizes('Inlier Concordance', 'ijiFig', ...
                this.mahalanobisConcordanceIdx, 'west++');
        end

        function rcScoresBySizes(this)
            if isnan(this.robustConcordanceIdx)
                msg('Robust concordance not computed this time');
                return;
            end
            this.scoresBySizes('Robust inlier concordance', 'rcFig', ...
                this.robustConcordanceIdx, 'north++');
        end

        function cScoresBySizes(this)
            this.scoresBySizes('Concordance', 'cFig', ...
                this.concordanceIdx, 'north west++');
        end

        function similarityScoresBySizes(this)
            this.scoresBySizes('Mass+Distance Similarity', 'similarityFig', ...
                this.similarityIdx, 'north east+');
        end

        function done=scoresBySizes(this, scoreTitle, figProperty, ...
                scoreIdx, where)
            if ishandle(this.(figProperty))
                figure(this.(figProperty));
                return;
            end
            scores=cell2mat(this.data(:, scoreIdx));
            scores(isnan(scores))=0;
            sizes=cell2mat(this.data(:, this.szIdx));
            trainers=strcmp('yes',this.data(:, this.groupIdx));
            scores=scores(trainers);
            sizes=sizes(trainers);
            [names, clrs]=LabelBasics.NamesAndColors(...
                this.qf.tIdPerRow(:,1), this.idMap);
            CellTypes=this.data(:, this.nameIdx);
            CellTypes=CellTypes(trainers);
            assert(isequal(names', CellTypes));
            [theTitle, done]=this.getClassificationTitle;
            if done
                this.(figProperty)=Plots.ScoreSizes([], scores, sizes,...
                    names, theTitle, scoreTitle, this.fig, where, clrs);
            end
        end

        function fig=confusionChart(this)
            try
                if ~Gui.IsVisible(this.confusionFig)
                    [~, this.confusionFig]=QfHiDM.ConfusionChart(this.qf);
                    n=this.fig.Name;
                    prfx=QfTable.FIG_NAME;
                    if startsWith(n, prfx)
                        n=['Confusion chart' n(length(prfx)+1:end)];
                    else
                        n='Confusion chart';
                    end
                    this.confusionFig.Name=n;
                    SuhWindow.Follow(this.confusionFig, ...
                        this.fig, 'south east+', true);
                    if Gui.IsVisible(this.fig)
                        SuhWindow.SetFigVisible(this.confusionFig);
                    end
                else
                    figure(this.confusionFig);
                end
                fig=this.confusionFig;
            catch ex
                ex.getReport
            end
        end

        function browseFalsePosNeg(this)
            if ~isempty(this.falsePosNegFig) ...
                    && ishandle(this.falsePosNegFig)
                figure(this.falsePosNegFig);
                return;
            end
            set(0, 'CurrentFigure', this.fig);
            busy=Gui.ShowBusy(this.fig, ...
                Gui.YellowSmall('Finding false positives & negatives'),...
                'wayneMoore1.png', .65, false);
            if isempty(this.externalArgs)
                descripts={'cell type', 'cells'};
            else
                try
                    descripts=this.externalArgs.class_descriptions;
                catch
                    descripts='';
                end
            end
            if isempty(busy)
                locatingFig=[];
            else
                locatingFig={this.fig, 'north east+', true};
            end
            [fig_,that]=FalsePositiveNegative.Plot([0 1], ...
                this.qf.falsePosNegFile, [], 2, descripts, ...
                'sample', false, false, true, locatingFig);
            if ~isempty(busy)
                figure(fig_);
            end
            that.predictions=this.getPredictionsOfThese;
            that.fcnMoreHtml=@()getFalsePosNegMatrixHtml(this);
            this.falsePosNegFig=fig_;
            Gui.AddSvgToToolBar(fig_);
            Gui.HideBusy(this.fig, busy, true);
        end
        
        
        function html=getFalsePosNegMatrixHtml(this)
            html=FalsePositiveNegative.MatrixHtml(this.qf);
        end

        function addToWebPage(this, fileName)
            visible=Gui.IsVisible(this.fig);
            this.doHistQF(visible);
            if this.qf.areEqual
                this.doHistF(visible);
            end
            this.browse('', false, false, fileName, true);
        end

        function saveSortedTable(this)
            T=this.sortTable.getSortedTableData(this.data);
            props=BasicMap.Global;
            prop1='QfTable.SaveFldr';
            prop2='QfTable.SaveFile';
            file=TableBasics.UiPutFile(prop1, prop2, props, ...
                'QFMatch_table.xls');
            if ~isempty(file) 
                if endsWith(lower(file), '.txt')
                    writetable(T, file, 'Delimiter', '\t');
                else
                    writetable(T, file);
                end
                File.OpenFolderWindow(file, 'QfTable.OpenFolder');
            end
        end

        function integrateClassification(this, showNow)
            busy=Gui.ShowBusy(this.fig, Gui.YellowH2(...
                'Integrating results with an Excel file'), ...
                'Genie.png', 1.1, true);
            try
                [T, this.qf.summaryFile, this.qf.classifierName, ...
                    this.qf.dataSetName, this.qf.testCaseName, ~, viewToo]...
                    =ClassificationTable.IntegrateResults(this, true, ...
                    '', '', ... %always give chance for NEW file
                    this.qf.classifierName,...
                    this.qf.dataSetName, this.qf.testCaseName);
                if ~isempty(T)
                    this.syncProperties;
                    this.summaryT=T;
                    this.saveClassification();
                    if showNow && viewToo
                        colorTable=table( ...
                            this.data(strcmp(this.data(:, this.groupIdx), ...
                            'yes'), this.nameIdx), ...
                            this.data(strcmp(this.data(:, this.groupIdx), ...
                            'yes'), this.symIdx));
                        try
                            %show MLP first
                            [comparisonTable, dataSetName]=...
                                ClassificationTable.ComparisonsTable( ...
                                T, this.qf.dataSetName, this.qf.testCaseName, ...
                                'MLP', colorTable);
                        catch
                            %MLP does NOT exist (I guess)
                            [comparisonTable, dataSetName]=...
                                ClassificationTable.ComparisonsTable( ...
                                T, this.qf.dataSetName, this.qf.testCaseName, ...
                                this.qf.classifierName, colorTable);
                        end
                        args=struct(...
                            'plotMdnMn', true, ...
                            'inlierConcordanceName', false, ...
                            'similarityName', false, ...
                            'genentechName', 'ESHGHI');
                        ClassificationTable.Html(comparisonTable, ...
                            dataSetName, true, this.fig, ...
                            'east+', this.qf.summaryFile, args);
                    end
                end
            catch ex
                Gui.MsgException(ex);
            end
            Gui.HideBusy(this.fig, busy, true);
        end
        
        function browse(this, matrixHtml, browseNow, ...
                addHeader, fileName, accumulate)
            if nargin<6
                accumulate=false;
                if nargin<5
                    fileName='';
                    if nargin<4
                        browseNow=true;
                        if nargin<3
                            addHeader=true;
                        end
                    end
                end
            end
            fh=this.fHistFig;
            fpn=this.falsePosNegFig;
            cc=this.confusionFig;
            try
                %earlier versions lack areEqual field in struct
                if ~this.qf.areEqual
                    fh=[];
                    fpn=[];
                    cc=[];
                end
            catch
            end
            predTable=[];
            try %not sure about every scenario
                %predTable=this.seePredictionOfThese; % =if ALWAYS want predictions (not sure)
                predTable=this.predictionsTable;
                if ~isempty(predTable)
                    predTable=predTable.sortTable.jtable;
                end
            catch ex
                ex.getReport
            end
            set(0, 'CurrentFigure', this.fig);
            QfTable.ToHtml(...
                this.sortTable.jtable, ...
                fh, ...
                this.qHistFig,...
                fpn,...
                predTable,...
                cc,...
                this.app, ...
                Html.remove(matrixHtml),...
                addHeader, ...
                browseNow, ...
                fileName, ...
                accumulate);
        end

        function changeName(this, priorName, newName, id, trainingSet, saveNow)
            if nargin<6
                saveNow=false;
                if nargin<5
                    trainingSet=false;
                end
            end
            t=this.sortTable.jtable;
            nameCol=t.convertColumnIndexToView(this.nameIdx-1);
            groupCol=t.convertColumnIndexToView(this.groupIdx-1);
            idCol=t.convertColumnIndexToView(this.idIdx-1);
            if isempty(trainingSet)
                matchGroup=-[];
            elseif trainingSet
                matchGroup='yes';
            else
                matchGroup='no';
            end
            nRows=t.getRowCount;
            changed=changeTable(priorName);
            if ~changed
                if contains(priorName, '_')
                    tryName=strrep(priorName, '_', '');
                    changed=changeTable(tryName);
                end
            end
            if changed
                nId=str2double(id);
                if isempty(trainingSet) || trainingSet
                    N=length(this.qf.tNames);
                    for i=1:N
                        if isequal(this.qf.tNames{i}, priorName)...
                                && this.qf.tIds(i)==nId
                            this.qf.tNames{i}=newName;
                            break;
                        end
                    end
                end
                if isempty(trainingSet) || ~trainingSet
                    N=length(this.qf.sNames);
                    for i=1:N
                        if isequal(this.qf.sNames{i}, priorName) ...
                                && this.qf.sIds(i)==nId
                            this.qf.sNames{i}=newName;
                            break;
                        end
                    end
                end
                if saveNow
                    this.save( this.qf, this.saveFile);
                end
            end
            
            function ok=changeTable(pName)
                ok=false;
                fprintf('Changing table %s to %s', pName, newName);
                for r=0:nRows-1
                    prev=t.getValueAt(r, nameCol);
                    if isequal(prev, pName)
                        id2=t.getValueAt(r, idCol);
                        if isequal(id, id2)
                            if isempty(trainingSet)
                                change=true;
                            else
                                group=t.getValueAt(r, groupCol);
                                change=isequal(matchGroup, group);
                            end
                            if change
                                t.setValueAt(newName, r, nameCol);
                                fprintf('CHANGED\n');
                                ok=true;
                                return;
                            end
                        end
                    end
                end
                fprintf('not CHANGED\n');
            end
        end
        
        function ok=doHistQF(this, visible)
            if nargin<2
                visible=true;
            end
            ok=true;
            if ishandle(this.qHistFig)
                if visible
                    figure(this.qHistFig);
                end
            else
                fig_=Gui.NewFigure(true, 'off');
                fig_.Color='white';
                this.qHistFig=fig_;
                Gui.Resize(fig_, this.histResize);
                ax=Gui.Axes(fig_);
                [qfData, avgMdn, avgMn]=this.getData;
                if isempty(qfData)
                    ok=false;
                    return;
                end                    
                histogram(ax, qfData, length(unique(qfData)));
                set(ax, 'units', 'normalized', 'Position', ...
                    [.1 .16 .8 .7]);
                Gui.StretchUpper(ax, @ylim, .1);
                Gui.StretchUpper(ax, @xlim, .1);
                xlabel(ax, '% mass+distance similarity', 'FontName', 'Arial')
                xlim(ax, [-5 105]);
                ylabel(ax, '# of subset matches', 'FontName', 'Arial')
                figName='% similarity ^{(mass+distance)}';
                this.addTitleToHist(fig_, ax, ...
                    figName, avgMdn, avgMn);
                this.addFalsePosNegButton(fig_);
                drawnow;
                if ~isempty(this.fHistFig) && ishandle(this.fHistFig)
                    where='west++';
                    followed=this.fHistFig;
                else                    
                	followed=this.fig;
                    where='south west++';
                end
                SuhWindow.Follow(this.qHistFig, followed, where, true);
                if visible
                    SuhWindow.SetFigVisible(fig_);
                end                
            end            
        end
        
        function addTitleToHist(this, fig, ax, score, avgMdn, avgMn)
            ttl=[score ', median/mean=\color{blue}' ...
                String.encodeRounded(avgMdn, 1) ...
                '\color{black}/\color{blue}'...
                String.encodeRounded(avgMn, 1) '\color{black}% '];
            set(fig, 'name', String.RemoveTex(score), 'NumberTitle', 'off');
            if ~isempty(this.matchesStr)
                title(ax, {ttl, this.matchesStr},  'FontName', 'Arial');
            elseif this.unmatched>0
                title(ax, {ttl, [ ...
                    num2str(this.R-this.unmatched) ' subsets matched, ' ...
                    '\color{red}' num2str(this.unmatched) ...
                    ' \color{black}NOT matched ']},'FontName', 'Arial');
            else
                title(ax, {ttl, [ ...
                    num2str(this.R) ' matches']}, ...
                    'FontName', 'Arial');
            end
        end
        
        function ok=doHistF(this, visible)
            ok=true;
            if nargin<3
                if nargin<2
                    visible=true;
                end
            end
            [fData, avgMdn, avgMn]=this.getData(0);
            if ishandle(this.fHistFig)
                if visible
                    figure(this.fHistFig);
                end
            else
                fig_=Gui.NewFigure(true, 'off');
                fig_.Color='white';
                this.fHistFig=fig_;
                Gui.Resize(fig_, this.histResize);
                ax=Gui.Axes(fig_);
                if isempty(fData)
                    ok=false;
                    return;
                end
                histogram(ax, fData, length(unique(fData)));
                xlim(ax, [-5 105]);
                ylabel(ax, '# of subset matches', 'FontName', 'Arial')
                set(ax, 'xlim', [0 100]);
                if this.qf.areEqual
                    this.addTitleToHist(fig_, ax, ...
                        '% Overlap ^{(F1-Score)}', avgMdn, avgMn);
                    xlabel(ax, '% overlap via F-measure', 'FontName', 'Arial')
                else
                    this.addTitleToHist(fig_, ax, ...
                        '% Overlap ^{(F-measure on probability bins)}', avgMdn, avgMn);
                    xlabel(ax, '% overlap via F-measure on probability bins', 'FontName', 'Arial')
                end
                set(ax, 'units', 'normalized', 'Position', ...
                    [.1 .16 .8 .7]);
                Gui.StretchUpper(ax, @ylim, .1);
                this.addFalsePosNegButton(fig_);
                if visible
                    drawnow;
                end
                if ~isempty(this.qHistFig) && ishandle(this.qHistFig)
                    where='east++';
                    followed=this.qHistFig;
                else                    
            	    followed=this.fig;
                    where='south++';
                end
                SuhWindow.Follow(this.fHistFig,followed, where, true);
                if visible
                    SuhWindow.SetFigVisible(fig_);
                end                
            end
            if visible
                if avgMdn<65
                    Gui.Laughter;
                elseif avgMdn >94
                    Gui.Applause;
                end
            end
        end
        
        function ok=doHistConcordance(this, visible)
            ok=true;
            if nargin<2
                visible=true;
            end
            if ishandle(this.cHistFig)
                if visible
                    figure(this.cHistFig);
                end
            else
                [fData, avgMdn, avgMn]=this.getData(2);
                if isempty(fData)
                    msgWarning(Html.WrapHr(['Robust inlier concordance'...
                        ' is<br>not calculated for this match...']));
                    ok=false;
                    return;
                end
                fig_=Gui.NewFigure(true, 'off');
                fig_.Color='white';
                this.cHistFig=fig_;
                Gui.Resize(fig_, this.histResize);
                ax=Gui.Axes(fig_);
                H=histogram(ax, fData, length(unique(fData)));
                xlabel(ax, '% concordance', 'FontName', 'Arial')
                xlim(ax, [-5 105]);
                ylabel(ax, '# of subset matches', 'FontName', 'Arial')
                set(ax, 'xlim', [0 100]);
                this.addTitleToHist(fig_, ax, ...
                    '% Concordance ^{(Jaccard index)}', avgMdn, avgMn);
                set(ax, 'units', 'normalized', 'Position', ...
                    [.1 .16 .8 .7]);
                Gui.StretchUpper(ax, @ylim, .1);
                this.addFalsePosNegButton(fig_);
                drawnow;
                if ~isempty(this.fHistFig) && ishandle(this.fHistFig)
                    where='east++';
                    followed=this.qHistFig;
                else                    
                	followed=this.fig;
                    where='south east++';
                end
                SuhWindow.Follow(this.cHistFig, followed, where, true);
                if visible
                    SuhWindow.SetFigVisible(fig_);
                end                
            end
        end
        
        function ok=doHistMahalanobisConcordance(this, visible)
             ok=true;
            if nargin<2
                visible=true;
            end
            if ishandle(this.mcHistFig)
                if visible
                    figure(this.mcHistFig);
                end
            else
                [fData, avgMdn, avgMn]=this.getData(3);
                if isempty(fData)
                    msgWarning(Html.WrapHr(['Robust inlier concordance'...
                        ' is<br>not calculated for this match...']));
                    ok=false;
                    return;
                end
                fig_=Gui.NewFigure(true, 'off');
                fig_.Color='white';
                this.mcHistFig=fig_;
                Gui.Resize(fig_, this.histResize);
                ax=Gui.Axes(fig_);
                H=histogram(ax, fData, length(unique(fData)));
                xlabel(ax, 'Central similarity', 'FontName', 'Arial')
                xlim(ax, [-5 105]);
                ylabel(ax, '# of subset matches', 'FontName', 'Arial')
                set(ax, 'xlim', [0 100]);
                this.addTitleToHist(fig_, ax, ...
                    '% Central similarity ^{(Inlier Concordance)}', avgMdn, avgMn);
                set(ax, 'units', 'normalized', 'Position', ...
                    [.1 .16 .8 .7]);
                Gui.StretchUpper(ax, @ylim, .1);
                this.addFalsePosNegButton(fig_);
                drawnow;
                if ~isempty(this.fHistFig) && ishandle(this.fHistFig)
                    where='east++';
                    followed=this.qHistFig;
                else                    
                	followed=this.fig;
                    where='south east++';
                end
                SuhWindow.Follow(this.mcHistFig, followed, where, true);
                if visible
                    SuhWindow.SetFigVisible(fig_);
                end                
            end
        end

        function ok=doHistRobustConcordance(this, visible)
            ok=true;
            if nargin<2
                visible=true;
            end
            if ishandle(this.rcHistFig)
                if visible
                    figure(this.rcHistFig);
                end
            else
                [fData, avgMdn, avgMn]=this.getData(4);
                if isempty(fData)
                    msgWarning(Html.WrapHr(['Robust inlier concordance'...
                        ' is<br>not calculated for this match...']));
                    ok=false;
                    return;
                end
                fig_=Gui.NewFigure(true, 'off');
                fig_.Color='white';
                this.rcHistFig=fig_;
                Gui.Resize(fig_, this.histResize);
                ax=Gui.Axes(fig_);
                H=histogram(ax, fData, length(unique(fData)));
                xlabel(ax, '% robust inlier concordance', 'FontName', 'Arial')
                xlim(ax, [-5 105]);
                ylabel(ax, '# of subset matches', 'FontName', 'Arial')
                set(ax, 'xlim', [0 100]);
                this.addTitleToHist(fig_, ax, ...
                    '% Robust Inlier Concordance ^{(Jaccard index)}', avgMdn, avgMn);
                set(ax, 'units', 'normalized', 'Position', ...
                    [.1 .16 .8 .7]);
                Gui.StretchUpper(ax, @ylim, .1);
                this.addFalsePosNegButton(fig_);
                drawnow;
                if ~isempty(this.fHistFig) && ishandle(this.fHistFig)
                    where='east++';
                    followed=this.qHistFig;
                else                    
                	followed=this.fig;
                    where='south east++';
                end
                SuhWindow.Follow(this.rcHistFig, followed, where, true);
                if visible
                    SuhWindow.SetFigVisible(fig_);
                end                
            end
        end
        
        function addFalsePosNegButton(this, fig_, lbl, width)
            if nargin<4
                width=.2;
                if nargin<3
                    lbl='See false +/-';
                end
            end
            if isempty(this.qf.falsePosNegs)
                return;
            end
            uicontrol(fig_, 'style', 'pushbutton','String', lbl,...
                'Units', 'normalized', ...
                'FontWeight', 'bold', ...
                'ForegroundColor', 'blue',...
                'BackgroundColor', [1 1 .80],...
                'Position',[1-(width+.01), .015, width, .071],...
                'Callback', @(btn,event) browseFalsePosNeg(this));
            uicontrol(fig_, 'style', 'pushbutton', ...
                        'String', 'Confusion chart',...
                        'Units', 'normalized', ...
                        'FontWeight', 'bold', ...
                        'ForegroundColor', 'blue',...
                        'BackgroundColor', [1 1 .80],...
                        'Position',[.03, .015, width*1.33, .071],...
                        'Callback', @(btn,event) confusionChart(this));
        end
        
        function QF=save(this, qf, file, warnNotSaved) 
            if isstruct(qf) || isempty(file)
                QF=qf;
                if nargin>3 && ~isempty(warnNotSaved)
                    msgWarning(['<html><b>This is a previous QFMatch' ...
                        '</b><hr><br>To permanently save ' warnNotSaved ...
                        '<br>re-run this QFMatch and redo ....<br>(sorry)'])
                end
                return;
            end
            QF.tIds=qf.tIds;
            QF.sIds=qf.sIds;
            QF.tNames=qf.tNames;
            QF.tSizes=qf.tSizes;
            QF.sSizes=qf.sSizes;
            QF.matrixHtml=qf.matrixHtml;
            QF.matrixHtml_browser=qf.matrixHtml_browser;
            QF.matrixHtmlOverlap=qf.matrixHtmlOverlap;
            QF.matrixHtmlOverlap_browser=qf.matrixHtmlOverlap_browser;
            QF.overlap_tNames=qf.overlap_tNames;
            QF.overlap_sNames=qf.overlap_sNames;
            QF.rowsI=qf.rowsI;
            QF.colsI=qf.colsI;
            QF.matrixOverlap=qf.matrixOverlap;
            QF.qfTableData=this.data;
            QF.unmatched=this.unmatched;
            QF.matchesStr=this.matchesStr;
            QF.tClrs=this.tClrs;
            QF.areEqual=qf.areEqual;
            if qf.areEqual
                qf.computeConfusion;
            end
            QF.confusionMatrix=qf.confusionMatrix;
            QF.confusionLabels=qf.confusionLabels;
            QF.falsePosNegs=qf.falsePosNegs;
            qf.relocateFalsePosNegFiles(file, true);
            QF.falsePosNegFile=qf.falsePosNegFile;
            QF.falsePosNegSubsetsFile=qf.falsePosNegSubsetsFile;
            QF.falsePosNegCnts=qf.falsePosNegCnts;
            QF.falsePosCulprits=qf.falsePosCulprits;
            QF.falseNegCulprits=qf.falseNegCulprits;
            QF.sNames=qf.sNames;
            QF.matches=qf.matches;
            QF.columnNames=qf.columnNames;
            QF.densityBars=qf.densityBars;
            QF.fMeasuringUnmerged=qf.fMeasuringUnmerged;
            [QF.tUnmatched, QF.tN, QF.sUnmatched, QF.sN]=qf.getMatchCounts;
            if ~isempty(qf.falsePosEvents)
                QF.predictions=SuhPredictions(qf);
                QF.predictions.clearMatchObject;
            else
                QF.predictions=[];
            end
            QF.tMdns=qf.tMdns;
            QF.tCompMdns=qf.tCompMdns;
            QF.sMdns=qf.sMdns;
            QF.sCompMdns=qf.sCompMdns;
            QF.classifierName=qf.classifierName;
            QF.dataSetName=qf.dataSetName;
            QF.testCaseName=qf.testCaseName;
            QF.summaryFile=qf.summaryFile;
            QF.saveFile=file;
            QF.matrixFinal=qf.matrixFinal;
            QF.sMergedIds=qf.sMergedIds;
            QF.tMergedIds=qf.tMergedIds;
            QF.fMeasuringMerged=qf.fMeasuringMerged;
            if nargin>2 && ~isempty(file)
                this.saveFile=file;
                QF.densityBars.app=[];%don't save the APP ... infinite loop
                save(file, 'QF');
                QF.densityBars.app=BasicMap.Global;
            end
        end

        function QF=saveClassification(this)
            saveFile_=this.saveFile;
            if isempty(this.saveFile)
                saveFile_=this.qf.saveFile;
            end
            if ~isempty(saveFile_)
                try
                    load(saveFile_, 'QF');
                    QF.classifierName=this.qf.classifierName;
                    QF.dataSetName=this.qf.dataSetName;
                    QF.testCaseName=this.qf.testCaseName;
                    QF.summaryFile=this.qf.summaryFile;
                    save(saveFile_, 'QF');
                catch ex
                    ex.getReport
                end
            end
        end

        function addSuffixToFigs(this, suffix)
            renameFig(this.fig);
            renameFig(this.qHistFig);
            renameFig(this.fHistFig);
            renameFig(this.cHistFig);
            renameFig(this.mcHistFig);
            renameFig(this.rcHistFig);
            
            function renameFig(fig)
                if ~isempty(fig)
                    set(fig, 'name', [ get(fig, 'name') ' ' suffix])
                end
            end
        end
        
        function closeFigs(this)
            Gui.CloseFig(this.fig);
            Gui.CloseFig(this.qHistFig);
            Gui.CloseFig(this.fHistFig);
            Gui.CloseFig(this.cHistFig);
            Gui.CloseFig(this.mcHistFig);
            Gui.CloseFig(this.rcHistFig);
            Gui.CloseFig(this.f1Fig);
            Gui.CloseFig(this.ijiFig);
            Gui.CloseFig(this.rcFig);
            Gui.CloseFig(this.cFig);
            Gui.CloseFig(this.similarityFig);
            
            Gui.CloseFig(this.falsePosNegFig);
        end

        function [sIds, tIds, sNames, tNames]=getUnmatched(this, str2Num)
            if nargin<2
                str2Num=true;
            end
            N=length(this.data);
            if str2Num
                tIds=[];
                sIds=[];
            else
                tIds={};
                sIds={};
            end
            tNames={};
            sNames={};
            for i=1:N
                if isnan(this.data{i, this.similarityIdx})
                    isT=strcmpi(this.data{i, this.groupIdx}, 'yes');
                    if isT
                        if str2Num
                            tIds(end+1)=str2double(this.data{i, this.idIdx});
                        else
                            tIds{end+1}=this.data{i, this.idIdx};
                        end
                        tNames{end+1}=this.data{i, this.nameIdx};
                    else
                        if str2Num
                            sIds(end+1)=str2double(this.data{i, this.idIdx});
                        else
                            sIds{end+1}=this.data{i, this.idIdx};
                        end
                        sNames{end+1}=this.data{i, this.nameIdx};
                    end
                end
            end
        end

        function [tUnmatched, tN, sUnmatched, sN, nMatchGroups]=getMatchCounts(this)
            try
            [tUnmatched, tN, sUnmatched, sN]=this.qf.getMatchCounts;
            catch
                tUnmatched=this.qf.tUnmatched;
                tN=this.qf.tN;
                sUnmatched=this.qf.sUnmatched;
                sN=this.qf.sN;    
            end
            nMatchGroups=length(this.qf.matches);
        end
        
        function [meanSimilarity, meanOverlap, missingSubsets, ...
                newSubsets, medianSimilarity, medianOverlap]...
                =getSummary(this)
            [~,medianSimilarity, meanSimilarity]=this.getData(1);
            [~,medianOverlap, meanOverlap]=this.getData(0);
            [missingSubsets, ~, newSubsets]=this.getMatchCounts; 
        end
        
        
        function rowIdxs=getRankSortedRowIndexs(this, rankings)
            if isempty(rankings)
                matchData=cell2mat(this.data(:, 3)); % column 3 is match
                rowIdxs=find(matchData==0);
                return;
            end
            rankData=cell2mat(this.data(:, 10)); % column 10 is rank
            l=ismember(rankData, rankings);
            
            [~,I]=sort(rankData(l));
            idxs=find(l);
            rowIdxs=idxs(I);
            %this.data(rowIdxs,10)'
        end
        
        function fields=getFields(this, rowIndex)
            fields=QfTable.RowToFields(this.data(rowIndex,:), this);
        end

        function associateGatingTree(this, gt)
            set(this.fig, 'UserData', this);
            this.gt=gt;
            gt.otherFigs{end+1}=this.fig;
        end        
    end
    
    methods(Static)
        function fields=RowToFields(row, this)
            if nargin<2
                fields.rank=row{10};
                fields.name=row{2};
                fields.similarity=row{4};
                fields.overlap=row{5};
                fields.trainingSet=row{6};
                fields.frequency=row{7};
                fields.count=row{8};
                fields.id=row{11};
                fields.matches=row{3};
                fields.row=row{1};
                fields.symbol=row{9};
            else
                fields.rank=row{this.rankIdx};
                fields.name=row{this.nameIdx};
                fields.similarity=row{this.similarityIdx};
                fields.overlap=row{this.overlapIdx};
                fields.trainingSet=row{this.groupIdx};
                fields.frequency=row{this.freqIdx};
                fields.count=row{this.szIdx};
                fields.id=row{this.idIdx};
                fields.matches=row{3};
                fields.row=row{1};
                fields.symbol=row{this.symIdx};
            end
        end

        function Save(qfTable, qf, file)
            qfTable.save(qf, file);
        end
        
        function QF=Load(file, complain, tData, tIdPerRow, sIdPerRow, sData)
            QF=[];
            try
                if exist(file, 'file')
                    load(file, 'QF');
                    if ~isempty(QF)
                        if ~isfield(QF, 'fMeasuringUnmerged')
                            QF.fMeasuringUnmerged=true;
                        end
                        if ~isfield(QF, 'classifierName')
                            QF.classifierName=[];
                        end
                        if ~isfield(QF, 'dataSetName')
                            QF.dataSetName=[];
                        end
                        if ~isfield(QF, 'testCaseName')
                            QF.testCaseName=[];
                        end
                        if ~isfield(QF, 'saveFile')
                            QF.saveFile=file;
                        end
                        if ~isfield(QF, 'densityBars')
                            QF.densityBars=[];
                        else
                            QF.densityBars.app=BasicMap.Global;
                        end
                        if ~isfield(QF, 'tMdns')
                            QF.tMdns=[];
                            QF.tCompMdns=[];
                            QF.sMdns=[];
                            QF.sCompMdns=[];
                        end
                        if nargin>2
                            QF.tData=tData;
                            if nargin>3
                                QF.tIdPerRow=tIdPerRow;
                            end
                            if nargin>4
                                if isequal(size(tIdPerRow), ...
                                        size(sIdPerRow))
                                    QF.sIdPerRow=sIdPerRow;
                                end
                            end
                            if ~isfield(QF, 'probability_bins')
                                QF.probability_bins=[];
                            end
                            if nargin>5
                                QF.sData=sData;
                            end
                        end
                    end
                end    
            catch ex
                if nargin==1 || complain
                    ex.getReport
                end
            end
        end
        

        function [data, names, widths, tips, unmatched, plotIdx, freqIdx,...
                rankIdx, symIdx, matchesStr, numIdx, nameIdx, matchesIdx, ...
                similarityIdx, overlapIdx, idIdx, szIdx, concIdx, ...
                mahalanobisConcIdx, robustConcIdx, idMap]=Contents(qf, clrs, pu)            
            numIdx=1;
            nameIdx=2;
            matchesIdx=3;
            similarityIdx=4;
            overlapIdx=5;
            plotIdx=6;
            freqIdx=7;
            szIdx=8;
            symIdx=9;
            rankIdx=10;
            concIdx=11;
            mahalanobisConcIdx=12;
            idMap=Map;
            robustConcIdx=nan;
            try
                [data, unmatched, matchesStr, ~, colors]...
                    =qf.getTableData(clrs);  
                hasRobustConcordance=size(data,2)>12;
            catch ex
                matchesStr='';
                if isstruct(qf)
                    data=qf.qfTableData;
                    unmatched=qf.unmatched;
                    if isfield(qf, 'matchesStr')
                        matchesStr=qf.matchesStr;
                    end
                else
                    data=[];
                    unmatched=0;
                    ex.getReport
                end
                hasRobustConcordance=size(data,2)>13;
            end
            needsConcordance=qf.areEqual;
            widths=[3 0; 19 nan; 5 0; 4 -1; 4 -1; ...
                7 nan; 3 -1;  6 0; 3 nan; 5 0];
            names={'#', ...
                'Subset (class) name', ...
                Html.WrapSm('Mat-<br>ches '), ...
                Html.WrapSm('Similarity, earth <br>mover''s distance'), ...
                Html.WrapSm('% Overlap<br>(F1-Score)'),...
                Html.WrapSm('Train-<br>ing set?'), ...
                Html.WrapSm('<br>Freq.'), ...
                Html.WrapSm('# of <br>events'), ...
                '<html><font size="6"><font  color="blue">&bull;</font></font></html>', ...
                Html.WrapSm('Rank')};
            tips={...
                'The subset''s item # ',...
                'The name of the subset including match name if applicable`',...
                'The number of subsets in this match set',...
                'Percentage of mass+distance similarity',...
                'Percentage of overlap of cells (or bins)',...
                'yes for training set, no for test set ',...
                ['The frequency of the subset within this' ...
                ' selection of cells/events'],...
                'The # of cells/events for the subset ',...
                'The subset''s frequency-sized symbol',...
                'The ranking of the match set'};
            if needsConcordance
                widths=[widths; 4 -1; 4 -1];
                names=[names, ...
                    Html.WrapSm('% Concordance<br>(Jaccard index)'),...
                    Html.WrapSm('Central similarity<br>(inlier concordance)')];
                tips=[tips, ...
                    'Percentage of Jaccard index between subsets', ...
                    ['Percentage of Jaccard index with outliers ' ...
                    'removed between subsets']];
                if hasRobustConcordance %includes concordance
                    robustConcIdx=13;
                    widths(end+1,:)=[4 -1];
                    names{end+1}=Html.WrapSm('% Robust Inlier<br>Concordance');
                    tips{end+1}=['Percentage of Jaccard index with ' ...
                        'outliers ROBUSTLY removed'];
                    idIdx=14;
                else
                    idIdx=13;
                end
            else
                idIdx=11;
            end
            
            % now add subset id column
            widths(end+1, :)=[6 nan];
            names{end+1}=Html.WrapSm('Subset<br>ID');
            tips{end+1}='Numeric identifier of subset (class)';

            tN=length(qf.tIds);
            sN=length(qf.sIds);
            if size(data,2) < idIdx % subset ID needs adding
                for i=1:tN
                    data{i, idIdx}=num2str(qf.tIds(i));
                end
                for i=1:sN
                    data{tN+i, idIdx}=num2str(qf.sIds(i));
                end
            end
            if exist('colors', 'var')
                for i=1:tN+sN
                    name=data{i, nameIdx};
                    id=data{i, idIdx};
                    idMap.set(id, name);
                    idMap.set([id '.clr'], colors(i,:))
                    idMap.set([id '.idx'], i)
                end
            else
                for i=1:tN+sN
                    name=data{i, nameIdx};
                    id=data{i, idIdx};
                    idMap.set(id, name);
                    sym=data{i,symIdx};
                    clrIdx=find(sym=='#',1);
                    if clrIdx<1
                        clr=Gui.HslColor(i,N);
                    else
                        clr=zeros(1,3);
                        clr(1)=hex2dec(sym(clrIdx+1:clrIdx+2));
                        clr(2)=hex2dec(sym(clrIdx+3:clrIdx+4));
                        clr(3)=hex2dec(sym(clrIdx+5:clrIdx+6));
                        clr=clr/255;
                    end
                    idMap.set([id '.clr'], clr)
                    idMap.set([id '.idx'], i)
                end
            end
            if exist('pu', 'var') && ~isempty(pu)
                if ~isempty(pu.pb)
                    pu.pb.setValue(N);
                end
            end
        end

        function qft=RunWithGt(pidSupervised, supervisors, fg, ...
                gid, fcs, fcsIdxs, data, file, embedding, visible, pu)
            qft=[];
            pFig=get(0, 'currentFig');
            gt=fg.gt;
            sameSample=isequal(gt.tp.getParentFileId(gid), ...
                    gt.tp.getParentFileId(pidSupervised));
            if exist(file, 'file') && ~QfTable.DEBUG
                if sameSample
                    [gtData, ~, ~, gtIds]=...
                        QfHiDM.GetRequiredData2(fcs, ...
                        fcsIdxs, fg.gt, gid, pu, visible);
                    qf=QfTable.Load(file,false,gtData,gtIds);
                else
                    qf=QfTable.Load(file);
                end
                if ~isempty(qf.columnNames)
                    qft=QfTable(qf, qf.tClrs, [], pFig, visible);
                    qft.fncSelect=@(qf, isTeacher, qfIdx)select(qf, isTeacher, qfIdx);
                    qft.doHistQF(visible);
                    if ~isempty(qft.qf.falsePosNegs)
                        qft.doHistF(visible);
                    end
                    [~, sLbls]=supervisors.getQfTrained(embedding);
                    return;
                end
            end
            if isempty(fcsIdxs)
                return;
            end
            newPu=nargin<10||isempty(pu);
            if newPu
                path=BasicMap.Path;
                pu=PopUp('Computing mass+distance similarity', 'north', ...
                    'Note...', false, true, fullfile(path, 'genieSearch.png'));
                pu.setTimeSpentTic(tic);
            else
                pu.setText2('Computing mass+distance similarity');
            end
            [gtData, ~, columnNames, gtIds]=QfHiDM.GetRequiredData2(fcs, fcsIdxs, ...
                gt, gid, pu, visible);
            gtData=fcs.ensureLogicle(gtData, fcsIdxs);
            %overlap collisions distort final resort
            gtIds=LabelBasics.RemoveOverlap(gtIds);
            if isempty(gtData)
                if ~isempty(pu)
                    pu.close;
                end
                return;
            end
            resizeSupervised=false; % get the predictions table!!
            if sameSample
                if ~isequal(size(data), size(gtData))
                    if sum(fg.clusterContext.pidSampleSelectedRows)...
                            ==size(data,1)
                        resizeSupervised=true;
                        data=gtData;
                    end
                end
            end
            if isequal(gtData, data)
                matchStrategy=2;
            else
                matchStrategy=1;
            end
            gids=MatBasics.GetUniqueIdsAndCntIfGt0(gtIds);
            N=length(gids);
            gtNames=cell(1,N);
            for i=1:N
                id=num2str(gids(i));
                gtNames{i}=gt.tp.getDescription(id);
            end
            [sNames, sLbls, clrs]=supervisors.getQfTrained(embedding);
            if resizeSupervised
                sLbls=MatBasics.ExpandIfPossible( ...
                        fg.clusterContext.pidSampleSelectedRows, sLbls);
            end
            qf=run_HiD_match(gtData, gtIds, data, sLbls, ...
                'trainingNames', gtNames, ...
                'testNames', sNames, ...
                'matchStrategy', matchStrategy, 'log10', true, 'pu', pu);
            [~,~,tIdxForS]=qf.getMatches;
            qf.setColumnNames(columnNames)
            gtClrs=zeros(N, 3);
            for i=1:N
                clrIdx=find(tIdxForS==i, 1, 'first');
                if clrIdx>0
                    gtClrs(i,:)=clrs(clrIdx,:);
                else
                    gid=gids(i);
                    gtClrs(i,:)=gt.highlighter.getColor(num2str(gid));
                end
            end
            if ~isempty(qf)
                qft=QfTable(qf, gtClrs, [], pFig, visible);
                qft.fncSelect=@(qf, isTeacher, qfIdx)select(qf, isTeacher, qfIdx);
                qft.doHistQF(visible);
                if matchStrategy==2
                    qft.doHistF(visible);
                end
                supervisors.annotateQfTable(qft, ...
                    matchStrategy, size(data,2));
                if ~isempty(file) && ~exist(file, 'file')
                    qft.save(qf, file);
                end
            end
            if newPu
                pu.close;
            end
            function select(qf, isTeacher, qfIdx)
                [name, lbl, tLbl, sLbl]=QfHiDM.GetIds2(qf, isTeacher, qfIdx);
                fprintf('teacher=%d, tId=%d, name="%s"\n',...
                    isTeacher, tLbl, name);
                if ~isempty(supervisors.btns)
                    supervisors.btns.get( find(supervisors.btnLbls==sLbl, 1)-1).doClick
                end
                if ~isempty(supervisors.args.parameter_names) && ...
                        supervisors.args.roi_table>=2
                    if isTeacher
                        d=gtData(gtIds==lbl,:);
                    else
                        d=data(sLbls==lbl,:);
                    end
                    try
                        needToMake=isempty(supervisors.roiTable) ...
                            || ~ishandle(supervisors.roiTable.table.table.fig);
                    catch
                        needToMake=true;
                    end
                    if needToMake
                        supervisors.roiTable=Kld.Table(d, ...
                            supervisors.args.parameter_names, ...
                            supervisors.args.roi_scales, qft.fig, name);
                        Gui.Locate(supervisors.roiTable.table.table.fig, ...
                            qft.fig, 'east++', true, true);
                    else
                        supervisors.roiTable.refresh(d, name);
                    end
                end
            end
        end
        function qft=Run(fcs, fcsIdxs, gt, gid1, gid2, file, visible, pu)
            pFig=get(0, 'currentFig');
            qft=[];
            if exist(file, 'file')
                qf=QfTable.Load(file);    
                qft=QfTable(qf, qf.tClrs, [], pFig, visible);
                qft.doHistQF(visible);
                return;
            end
            if isempty(fcsIdxs)
                return;
            end
            newPu=nargin<8|| isempty(pu);
            if newPu
                path=BasicMap.Path;
                pu=PopUp('Computing mass+distance similarity', 'north', ...
                    'Note...', false, true, fullfile(path, 'genieSearch.png'));
                pu.setTimeSpentTic(tic);
            else
                pu.setText2('Computing mass+distance similarity table');
            end
            [data1, ids1, names1, clrs1]=go(gid1);
            if isempty(data1)
                return;
            end
            [data2, ids2, names2]=go(gid2);
            if isempty(data2)
                return;
            end
            matchStrategy=1;
            qf=run_HiD_match(data1, ids1, data2, ids2, ...
                'trainingNames', names1, 'testNames', names2, ...
                'matchStrategy', matchStrategy, 'log10', true, ...
                'pu', pu);
            if ~isempty(qf)
                qft=QfTable(qf, clrs1, [], pFig, visible);
                qft.doHistQF(visible);
                if ~isempty(file) && ~exist(file, 'file')
                    qft.save(qf, file);
                end
            end
            if newPu
                pu.close;
            end
            
            function [gtData, gtIds, names, clrs]=go(gid)
                [gtData, ~, ~, gtIds]=QfHiDM.GetRequiredData2(fcs, ...
                    fcsIdxs, gt, gid, pu, visible);
                if isempty(gtData)
                    return;
                end
                gids=MatBasics.GetUniqueIdsAndCntIfGt0(gtIds);
                hghl=gt.highlighter;
                N=length(gids);
                clrs=zeros(N, 3);
                names=cell(1,N);
                for i=1:N
                    id=num2str(gids(i));
                    names{i}=gt.tp.getDescription(id);
                    clrs(i,:)=hghl.getColor(id, Gui.HslColor(i, N));
                end
            end
        end        
        
        function [selectedLbls, selectedBtns]=GetSelectedLabels(btns, btnLbls)
            N=btns.size;
            selectedLbls=[]; selectedLbls(N) = 0;
            selectedBtns=cell(1,N);
            nSelected=0;

            for i=1:N
                btn=btns.get(i-1);
                if btn.isSelected
                    nSelected = nSelected+1;
                    selectedLbls(nSelected)=btnLbls(i);
                    selectedBtns{nSelected}=btn;
                end
            end

        selectedLbls = selectedLbls(1:nSelected);
        selectedBtns = selectedBtns(1:nSelected);
        end
        
        
         function [html, fileName, accumulate]...
                =ToHtml(jtable, fHist, qfHist,...
                falsePos, predTable, confusion,...
                props, additionalFooter,...
                addHeader, browseNow, fileName,...
                accumulate)
             if nargin<12
                 accumulate=false;
                 if nargin<11
                     fileName='';
                     if nargin<10
                         browseNow=true;
                         if nargin<9
                             addHeader=true;
                             if nargin<8
                                 additionalFooter=[];
                                 if nargin<7
                                     props=BasicMap.Global;
                                 end
                             end
                         end
                     end
                 end
             end
            prop='QfTable.Html.File';
            prop2='QfTable.Html.Folder';
            html='';
            figHtml='';
            try
                tableTitle=char(Gui.WindowAncestor(jtable).getTitle);
                mainFigIsVisible=true;
            catch 
                tableTitle=get(gcf, 'Name');
                mainFigIsVisible=false;
            end
            if addHeader
                btn=Gui.NewBtn(Html.WrapSmallBold(...
                    'Specify an<br>HTML file'), ...
                    @(h,e)putFile(), 'Put results into a specific file', ...
                    'file_open.png');
                cb=Gui.CheckBox(Html.WrapSmallBold(...
                    'Add if<br>file preexists'),  ...
                    false, props, 'QfMatch.AddHtml', [], ...
                    'Accumulate results of many matches in one html file');
                prop3='QfTable.Html.Title';
                props.set(prop3, tableTitle);
                [h1Header, cancelled]=inputDlg(...
                    struct('property', prop3,...
                    'properties', props, 'msg', ...
                    Gui.BorderPanel([], 2, 3, ...
                    'South', ['<html><hr><br>Enter a title for '...
                    'these match results...</html>'], ...
                    'West', btn, 'East', cb),...
                    'where', 'west'), 'Seeing match results in browser...', ...
                    'Describe this Hi-D match results');
                if cancelled
                    return;
                end
                accumulate=cb.isSelected;
                h1Header=['<hr><h2>' h1Header '</h2>'];
            else
                h1Header=['<hr><h2>' tableTitle '</h2>'];
            end
            if browseNow
                pu=PopUp('Generating webpage...');
            end
            html=SortTable.ToHtml(jtable);
            if ~isempty(predTable)
                if ~mainFigIsVisible || ~isempty(Gui.WindowAncestor(predTable))
                    %is table showing?
                    html=[html '<hr>' SortTable.ToHtml(predTable)];
                end
            end
            if isempty(fileName)
                fldr=tempdir;
            else
                fileName=File.ExpandHomeSymbol(fileName);
                fldr=fileparts(fileName);
                if isempty(fldr)
                    fldr=tempdir;
                end
            end
            if ~exist(fldr, 'dir')
                File.mkDir(fldr);
            end
            if ~isempty(fHist) ...
                    || ~isempty(qfHist) ...
                    || ~isempty(falsePos)
                histCount=0;
                figHtml=[figHtml ...
                    '<table border="0"><tr><td>'];
                figFile(fHist, 'F-Measure', .65);
                figFile(qfHist, 'QFMatch', .65);
                figHtml=[figHtml '</td></tr></table>'];
                if isempty(falsePos)
                    html=['<table border="1"><tr><td '...
                        'colspan="2" align="center">' ...
                        h1Header '</td></tr>'...
                        '<tr><td valign="top">'...
                        figHtml '</td><td valign="top" '...
                        'rowspan="2">' html '</td></tr>'...
                        '<tr><td valign="bottom">'...
                        '<center><small>(repeat of above '...
                        'histograms<br>for scrolling conv'...
                        'enience)</small></center>' figHtml ...
                        '</td></tr></table>' additionalFooter];
                else
                    img1=figHtml;
                    figHtml='';
                    figFile(falsePos, 'False +/-', .5);
                    html=['<table border="1"><tr><td '...
                        'colspan="2" align="center">' ...
                        h1Header '</td></tr>'...
                        '<tr><td valign="top">'...
                        img1 '</td><td valign="top" '...
                        'rowspan="2">' html '</td></tr>'...
                        '<tr><td valign="bottom">'...
                        '<center><small>False positive and'...
                        '<br>false negative)</small></center>' figHtml ...
                        '</td></tr></table>' ];
                end
            else
                html=[h1Header html ];
            end
            if ~isempty(confusion)
                figHtml='';
                figFile(confusion, 'confusion',.97);
                html=[html '' figHtml];
            end
            if ~isempty(additionalFooter)
                html=[html additionalFooter];
            end
            img1=Html.Img('epp.png', [], .10, true);
            img2=Html.Img('epp.png', [], .10, false);
            html=strrep(html, img2, img1);
            if isempty(fileName)
                if accumulate
                    fileName=props.get(prop);
                end
                if isempty(fileName)
                    fileName=File.SaveTempHtml(html);
                end
            end
            if accumulate
                if File.ExistsFile(fileName)
                    priorHtml=File.ReadTextFile(fileName);
                    html=[Html.remove(priorHtml) ...
                        html];
                end
            end
            File.SaveTextFile(fileName, ...
                ['<html>' html '</html>']);
            props.set(prop, fileName);
            props.set(prop2, fileparts(fileName));
            if browseNow
                web( fileName, '-browser');
                pu.close;
            end


            function putFile()
                defaultName=File.Time('QFMatch', 'html');
                [fldr2, file2]=uiPutFile(File.Documents, ...
                    defaultName, props, prop2,...
                    'Specify an HTML file....');
                if ~isempty(fldr2)
                    fileName=fullfile(fldr2, file2);
                end
            end
            
            function figFile(fig, filePrefix, scale)
                if ~isempty(fig) && (~mainFigIsVisible ...
                        || Gui.IsVisible(fig))
                    if isempty(fileName)
                        f='';
                    else
                        [~, f]=fileparts(fileName);
                        f=[f '_'];
                    end
                    histCount=histCount+1;
                    if histCount==2
                        br='</td></tr><tr><td>';
                    else
                        br='';
                    end
                    file=File.Time([f filePrefix], 'png');
                    saveas(fig, fullfile(fldr, file));
                    figHtml=[figHtml br Html.ImgXy(...
                        file, fldr, scale, true)];
                end
            end
        end
    end
end
