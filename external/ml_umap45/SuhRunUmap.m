classdef SuhRunUmap < handle
%SuhRunUmap aims to make the old run_umap more readable. It became
% unreadable between May 2019 and October 2021 due to rapid addition of
% easy end-user features and arguments ... a typical feature creep.
%
% The method for achieving readability is subsetting related code into
% function enclosed functions that access local variables of the top
% constructor function. Most of the data state in running UMAP does not
% need persistent class variables except for the original output arguments
% of run_umap (reduction, umap, clusterIdentifiers and extras).

%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause
    

    properties(Constant)
        CLUSTER_BORDER_CLR=[.15 .0 .3];
    end
    properties(SetAccess=private)
        reduced_data;
        umap;
        extras
        clusterIdentifiers;
        fncSaveAllClusters;
        fig;
        sprv;
    end
    
    methods(Static)
        function [reduction, umap, clusterIdentifiers, extras]=Go(varargin)
            this=SuhRunUmap(nargout, varargin{:});
            reduction=this.reduced_data;
            umap=this.umap;
            clusterIdentifiers=this.clusterIdentifiers;
            extras=this.extras;            
        end
    end
    
    methods 
        function saveOriginalClusters(this, unmatchedOnly)
            if nargin<1
                unmatchedOnly=false;
            end
            if ~isempty(this.fncSaveAllClusters)
                if unmatchedOnly
                    if isempty(this.extras.qfd)
                        warning('No QFMatch done!');
                    else
                        sIds=this.extras.qfd{1}.getUnmatched(true);
                        if isempty(sIds)
                            warning('QFMatch matched ALL clusters!');
                        else
                            feval(this.fncSaveAllClusters, sIds);
                        end
                    end
                else
                    feval(this.fncSaveAllClusters);
                end
            end
        end

        function this=SuhRunUmap(nArgOut, varargin)
            this.reduced_data=[];
            this.umap=[];
            this.extras=[];
            this.clusterIdentifiers=[];
            curPath=fileparts(mfilename('fullpath'));
            if ~CheckUmapFolder(curPath, 'FileBasics.m', true)
                return;
            end
            if ~CheckUmapFolder(curPath, UmapUtil.LocateMex)
                return;
            end
            try
                if length(varargin)==2 % looks odd
                    if startsWith('job_folder', lower(varargin{1}))
                        UmapUtil.RunJobs(varargin{2});
                        return;
                    end
                end
                [args, argued]=UmapUtil.Initialize(varargin{:});
            catch ex
                Gui.MsgException(ex);
                error(ex.message);
                return;
            end
            
            if ~isempty(args.job_folder)
                UmapUtil.RunJobs(args.job_folder);
            end
            if argued.contains('randomize')
                if ~args.randomize
                    warning(['randomize==true is faster preserving local/global structure'...
                        ' & ONLY changing south/north/south/east orientation']);
                end
            end
            args=UmapUtil.CheckArgs(args, argued);
            globals=BasicMap.Global;
            globals.currentJavaWindow=[];%in case odd prior exception
            if isempty(globals.propertyFile)
                homeFolder=globals.appFolder;
                try
                    props=fullfile(homeFolder, BasicMap.FILE);
                    globals.load(props);
                catch
                end
            end
            args.reduction=floor(etime(clock, datevec('Jan 16 2021','mmmm dd yyyy')));
            reduction=[]; mins=[]; maxs=[];
            sprvLabels={};
            csv_file_or_data=args.csv_file_or_data;
            extras=UMAP_extra_results;
            beQuiet=strcmpi(args.verbose, 'none');
            if ~beQuiet
                disp(UMAP.DISCLOSURE);
            end
            firstQf=[];
            MatBasics.WarningsOff;
            
            beGraphic=strcmpi(args.verbose, 'graphic');
            curAxes=[];
            beGraphic = false;
            
            %{
            if beGraphic
                %init top level variables set by prepareGui
                firstPlot=true;xLabel=[];yLabel=[];zLabel=[];
                askedToSeeOutputFolder=false;            
                try
                    [fig, curAxes, tb, btnClusterRoi, cbSyncKld, comboClu, ...
                        lblClu, btnNameRoi, btnUst]=prepareGui;
                catch ex
                    halt=~askYesOrNo(Html.WrapHr(...
                        ['If you are running the <br>"<b>New ' ...
                        'Desktop for MATLAB</b>"' ...
                        '<br>add-on then you MUST disable' ...
                        '<br> or uninstall it!<br><br>...' ...
                        '<hr>Continue <i>without SEEING plot?</i>']), ...
                        ' GUI environment problems..', [], false);
                    Gui.MsgException(ex);
                    if halt
                        throw(ex);
                    else
                        beGraphic=false;
                    end
                end
                if ~handleNoData
                    return;
                end
            end

            %}
            
            %init top level variables set by prepareData
            hClusterLabels=[];
            cntClusterLabels=0;
            jdHeatMap=[];
            hClusterBorders={};
            hMeasurement=[];
            hMeasurementIdx=0;
            hMeasurementLabel=[];
            
            inData=[];nRows=[];parameter_names=[]; nCols=[]; 
            labels=[]; labelCols=[]; testSetLabels=[]; matchedLbls=[];
            newSubsetIdxs=[];umap=[];sprv=[];sCols=[];
            inputDataIsCsv=ischar(args.csv_file_or_data) && ...
                endsWith(lower(args.csv_file_or_data), '.csv');
            if (ischar(args.csv_file_or_data) && ...
                    endsWith(lower(args.csv_file_or_data), '.wsp')) ...
                    || (iscell(args.csv_file_or_data) && ...
                    endsWith(lower(args.csv_file_or_data{1}), '.wsp'))
                try
                    if beGraphic
                        figFjb=fig;
                        set(get(fig, 'currentAxes'), 'XTick', [], ...
                            'YTick', []);
                        drawnow;
                        fjwBusy=title(Gui.Axes(fig),{...
                            'Opening FlowJo 10.x workspace',...
                            'for data and gating classifications'});
                    else
                        figFjb=[];
                    end
                    if args.n_components ~= 2
                        msgWarning(Html.SprintfHr(['Changing ' ...
                            'n_components from %d to 2 which <br>is' ...
                            'required by FlowJoBridge. '], args.n_components))
                        args.n_components=2;
                    end
                    [args.csv_file_or_data, args.parameter_names, ...
                        args.label_file, args.flowjo_tree,...
                        args.sample_offsets, args.conclude, ...
                        args.save_roi]...
                        =FlowJoTree.Read(...
                        args.csv_file_or_data,...
                        args.flowjo_columns, ...
                        args.flowjo_visible, ...
                        args.label_column==0, ...
                        false, figFjb, args.conclude, ...
                        'UMAP', args.flowjo_ask);
                    if ~isempty(args.label_file) ...
                            && exist(args.label_file, 'file')
                        args.buildLabelMap=false;
                    end
                    if isempty(args.csv_file_or_data)
                        if ~isempty(args.flowjo_tree) && args.flowjo_tree.isTreeVisible
                            msg(['<html>Select a gate from the <br>'...
                                'tree then click "<font color="blue">' ...
                                'Run UMAP</font>"<br>' ...
                                'at the bottom right.'], 8, 'east+');
                            args.flowjo_tree.addPipelineBtn(varargin);
                        else
                            MatBasics.RunLater(...
                                @(h,e)msg(['<html>No reduction will '...
                                'be done...</html>']), 1);
                        end
                        if beGraphic
                            close(fig);
                        end
                        return;
                    end
                    if isnan(args.std_outliers)
                        args.std_outliers=3;
                    end
                    if beGraphic
                        delete(fjwBusy);
                    end
                catch ex
                    Gui.MsgException(ex);
                    if beGraphic
                        close(fig);
                    end
                    
                    return;
                end
            end
            [ok, template_file]=prepareData;
            if ~ok
                return;
            end
            if islogical(args.match_webpage_file) ...
                    && args.match_webpage_file
                if ~inputDataIsCsv
                    warning(['If match_webpage_file==true THEN ' ...
                        'csv_file_or_data must be a file name']);
                else
                    [hWait, f]=fileparts(args.csv_file_or_data);
                    args.match_webpage_file=fullfile(hWait, [f '_webpages'],...
                        [f '.html']);
                    if ~argued.containsStartsWithI('false_positive_negative_plot')
                        args.false_positive_negative_plot=true;
                    end
                    if ~argued.containsStartsWithI('confusion_chart')
                        args.confusion_chart=true;
                    end
                end
            end
            if ischar(args.match_webpage_file) ...
                    && ~isempty(args.match_webpage_file)
                if args.match_webpage_reset
                    try
                        delete(args.match_webpage_file, 'file');
                    catch
                    end
                end
            end
            
            %init top level local variables set by prepareUmapSettings
            progressMatchType=0;method=[]; isSupervising=false;
            if ~prepareUmapSettings
                return;
            end
            
            %init top level local variables set by prepareLabels
            labelMap=[];hasLabels=false;nLabels=0;
            prepareLabels;
            
            %init top level local variables set by handleBadData
            if ~handleBadData
                return;
            end
            if beGraphic
                Gui.Trumpet;
            end

            %init top level function variables set by prepareAnnotations
            paramAnnotation=[];runAnnotation=[]; 
            pu=[];
            if beGraphic
                if ~ishandle(curAxes)
                    return;
                end
                simplicialIdx=1;
                simplicialSize=globals.adjustHighDef(.152,.05);
                busy=[];
                prepareAnnotations;
            end
            
            tick=tic;                
           
            %READY TO REDUCE Hi-D to Lo-D!!
            %init top level function variables set by reduceHiD2LoD
            reduction=[]; strReduction=[];pythonTemplate=[];...
                probability_bins=[];clusterIdentifiers=[];densityBars=[];
            
            reduceHiD2LoD;
            
            if ~isempty(reduction)
                %prepare top level function variables set by
                %saveShowReduction
                roiMap=[];legendRois=[];dbm=[];timeAnnotation=[];
                if beGraphic
                    Gui.BlowWindow(tb.jToolbar);
                end
                showAndSaveReduction;
            else
                if ~args.mlp_only
                    msg(Html.WrapHr(['Parameter reduction was cancelled ...'...
                        '<br>' Html.WrapBoldSmall(' (or not done)') ]));
                end
                if exist('pu', 'var') && isa(pu, 'PopUp')
                    pu.stop;
                    pu.dlg.dispose;
                end
            end
            
            if ~beGraphic
                if nArgOut<4 %don't delete figures if expecting them in extras
                    extras.closeMatchFigs;
                    extras.closeTreeFigs;
                end
            end
            globals.save;
            if beGraphic && ~isempty(reduction)
                roi=[];roiTable=[]; drawingRoi=false;
                seekingDataIsland=false;wbUp=true;
                lastRoiPos=[];dbmROI=[];jdDbm=[];                   
                finishPlot;
            end
            this.reduced_data=reduction;
            if ~isempty(umap)
                umap.args=args;
            end
            this.umap=umap;
            this.extras=extras;
            this.clusterIdentifiers=clusterIdentifiers;
            listenToPredictions;
            sampleRows=[];
            highlighted=[];
            if beGraphic && ishandle(fig)
                if ~isempty(args.match_webpage_file)
                    Gui.AskYesOrNoNonModal(...
                        ['<html>Browse webpage now?' ...
                        Html.FileTree(args.match_webpage_file) '</html>'],...
                        @browse, 9, 'south++', fig);
                end
                if ~isempty(args.highlighter_registry)
                    sampleRows=args.sample_rows;
                    sampleOffsets=args.sample_offsets;
                    if ~isempty(sampleOffsets)
                        sampleRowsFalse=false(1, length(sampleRows));
                    end
                    highlighted=Map;
                    feval(args.highlighter_registry, ...
                        @(g,o,r)hearHighlighting(g,o, r));
                end
            end
            if ~isempty(args.conclude)
                try
                    if beGraphic
                        ff=fig;
                        args.fig=ff;
                    else
                        ff=[];
                    end
                    if beGraphic
                        tb.setEnabled(false);
                    end
                    feval(args.conclude, reduction, sprv, args, ff);
                catch ex
                    ex.getReport
                end
            end
            if beGraphic && ~isempty(busy)
                tb.setEnabled(true);
                Gui.HideBusy(fig, busy);
                busy=[];
            end
            this.fncSaveAllClusters=@convertAllClusters;
            this.sprv=sprv;
            if beGraphic
                this.fig=fig;
            end

            function ok=hearHighlighting(gate, on, rows)
                if ~ishandle(curAxes)
                    ok=false;
                    return
                end
                ok=true;
                g=gate.gater;
                if nargin<3
                    rows=g.getSampleRows(gate);
                end
                wasHeld=ishold(curAxes);
                if ~wasHeld
                    hold(curAxes, 'on');
                end
                %fprintf('%s (%s) highlighting==%d\n', gate.name, gate.id, on);
                if highlighted.containsKey(gate.id)
                    H=highlighted.get(gate.id);
                    delete(H);
                end
                if isempty(sampleOffsets)
                    rows_=sampleRows&rows;
                else
                    rows_=sampleRowsFalse;
                    offsets=sampleOffsets.get(gate.id);
                    pg=gate;
                    while isempty(offsets)
                        pg=pg.getParent;
                        if pg.isSample
                            break;
                        end
                        offsets=sampleOffsets.get(pg.id);
                    end
                    if isempty(offsets)
                        ok=false;
                        return;
                    end
                    ending=offsets.sampleOffset+offsets.sampleSize-1;
                    rows_(offsets.sampleOffset:ending)=rows;
                end
                r=rows_(sampleRows);
                nHighlights=sum(r);
                Gui.RemoveTitleFlashStr(curAxes, gate.name);
                if nHighlights>0 && on
                    H=Gui.Flash(curAxes, reduction(r,:), ...
                        gate.name, false, 6, gate.highlightColor);
                    set(H, 'PickableParts', 'none');
                    highlighted.set(gate.id, H);
                end
                if ~wasHeld
                    hold(curAxes, 'off');
                end
                %N=length(g.highlightedGates)
            end
            
            function ok=browse(~, answ)
                ok=true;
                if contains(lower(answ), 'yes')
                    web(args.match_webpage_file, '-browser');
                end
            end

            function listenToPredictions
                if beGraphic && ~isempty(this.extras)
                    qft=this.extras.getMatchTable(5);
                    if ~isempty(qft) && ~isempty(qft.predictions)
                        qft.predictions.setSelectionListener(...
                            @notifyPredictions);
                    end
                    qft=this.extras.getMatchTable(3);
                    if ~isempty(qft)
                        qft.setPredictionListener(@notifyPredictions);
                    end
                    qft=this.extras.getMatchTable(4);
                    if ~isempty(qft)
                        qft.setPredictionListener(@notifyPredictions);
                    end
                end
            end
            
            function notifyPredictions(predictions)
                if ~isempty(predictions.selectedIds)
                    Gui.Flash(curAxes, ...
                        reduction(predictions.selectedData,:), ...
                        predictions);
                end
            end
            
            function finishPlot
                %handle region of interest stuff
                if exist('tb', 'var') && ishandle(fig)
                    set(fig,'WindowButtonDownFcn', @(h,e)wbd);
                    set(fig, 'WindowButtonUpFcn', @(h,e)wbu);
                    if ~isempty(tb)
                        tb.setEnabled(true);
                    end
                end
                setClusterDetail;
                MatBasics.RunLater(@(h,e)hideAnnotations, 7);
                if ~argued.contains('plot_title')
                    if inputDataIsCsv
                        [~,csvFile]=fileparts(csv_file_or_data);
                        if args.fast_approximation
                            word='\itapproximated\rm';
                        else
                            word='reduced';
                        end
                        ttlCsvFile=String.ToLaTex(csvFile);
                        if ~args.python
                            title(curAxes, [args.reduction_algorithm ...
                                ' ' word ...
                                ' \color{blue}' ttlCsvFile]);
                        else
                            title(curAxes, {...
                                [args.reduction_algorithm ...
                                ' ' word ' \color{blue}' ...
                                ttlCsvFile],...
                                '\color{magenta}(using Python)'});
                        end
                    end
                elseif ~isempty(args.plot_title)
                    title(curAxes, args.plot_title);
                end
            end
            
            function showAndSaveReduction
                if args.save_output
                    if isempty(args.output_folder)
                        if inputDataIsCsv
                            args.output_folder=fileparts(csv_file_or_data);
                        else
                            args.output_folder=fullfile(File.Documents, 'run_umap');
                        end
                        if beQuiet
                            warning('No output_folder given so using %s', ...
                                args.output_folder);
                        else
                            args.output_folder=File.GetDir(...
                                fullfile(File.Documents, 'run_epp'), ...
                                'run_umap.output_folder',...
                                'folder for run_umap output');
                            if isempty(args.output_folder)
                                args.save_output=false;
                            end
                        end
                    end
                    if args.save_output
                        outputFile=['umap_' args.reductionType...
                            args.output_suffix ];
                        csvFile=fullfile(args.output_folder, [outputFile '.csv']);
                        File.mkDir(args.output_folder);
                        File.WriteCsvFile(csvFile, reduction, {'UMAP_X', 'UMAP_Y'});
                    end
                end
                extras.timing=toc(tick);
                if beGraphic && ishandle(fig)
                    figure(fig);
                    if ~exist('pu', 'var')
                        pu=PopUp('Updating plot', 'north east', [], false);
                    end
                    delete(runAnnotation);
                    if strcmpi(method, 'Java') || strcmpi(method, 'C++') || strcmpi(method, 'MEX')
                        pu.pb.setString('All done');
                    end
                    updatePlot(reduction, true)
                    if args.save_output
                        set(paramAnnotation, 'visible','off')
                        drawnow;
                        pngFile=fullfile(args.output_folder, [outputFile '.png']);
                        set(paramAnnotation, 'visible','on')
                        Gui.SavePng(fig, pngFile);
                    end
                    if args.fast_approximation
                        timeAnnotation=Gui.TextBox(['Approximation time=\color{blue}' ...
                            String.MinutesSeconds(extras.timing)], fig, ...
                            'position', [.6 .01 .39 .05]);
                    else
                        timeAnnotation=Gui.TextBox(['Reduction time=\color{blue}' ...
                            String.MinutesSeconds(extras.timing)], fig, ...
                            'position', [.65 .01 .33 .05]);
                    end
                    if isempty(template_file) && args.ask_to_save_template
                        saveTemplate();
                    end
                end
                testBasicReduction;
                if ~beQuiet
                    fprintf('%s reduction finished (cost %s)\n', ...
                        args.reduction_algorithm, ...
                        String.MinutesSeconds(toc(tick)));
                end
                if args.fast_approximation
                    reportProgress(['Finished (fast approximation) ' strReduction]);
                else
                    reportProgress(['Finished ' strReduction]);
                end
                if ~beGraphic
                    if isSupervising
                        if ~beQuiet
                            pu=PopUp('Matching results', 'west+', [], false, [], [],...
                                false, args.parent_popUp);
                        else
                            pu=[];
                        end
                        if nArgOut>3
                            if ~beQuiet
                                disp('Setting supervisor labels');
                            end
                            nSupervisors=sprv.computeAndMatchClusters( ...
                                reduction, args.match_supervisors(1), pu);
                            if nSupervisors>0
                                extras.supervisorMatchedLabels=sprv.supervise(...
                                    reduction, false, args.match_supervisors(1));
                            end
                        end
                        doQfs(reduction);
                        if ~beQuiet
                            pu.close;
                        end
                    end
                end
                N=length(extras.qfd);
                for j=1:N
                    qfFig=extras.qfd{j}.fig;
                    set(qfFig, 'name', [ get(qfFig, 'name') ' ('...
                        String.encodeInteger(extras.qfd{j}.qf.matchTiming)...
                        ' secs)'])
                end
                if args.save_output
                    for j=1:N
                        qfd=extras.qfd{j};
                        outputFile=['similarity_histogram_' num2str(j) args.output_suffix];
                        pngFile=fullfile(args.output_folder, [outputFile '.png']);
                        Gui.SavePng(qfd.qHistFig, pngFile);
                        outputFile=['overlap_histogram_' num2str(j) args.output_suffix];
                        pngFile=fullfile(args.output_folder, [outputFile '.png']);
                        Gui.SavePng(qfd.fHistFig, pngFile);
                    end
                    N=length(extras.qft);
                    for j=1:N
                        outputFile=['qf_tree_' num2str(j) args.output_suffix];
                        pngFile=fullfile(args.output_folder, [outputFile '.png']);
                        Gui.SavePng(extras.qft(1).fig, pngFile);
                    end
                end
                if beGraphic
                    doLegendRoiButton;
                    if isempty(legendRois)
                        roiMap=Map;
                    else
                        roiMap=legendRois.map;
                    end
                end
                if (nArgOut>1 || ~isempty(args.save_template_file)) 
                    if beGraphic
                        [stuff, oldArgs]=umap.prepareForTemplate(curAxes, true);
                    else
                        [stuff, oldArgs]=umap.prepareForTemplate([], true);
                    end
                    if args.python
                        if  ~isempty(args.save_template_file)
                            [f1, f2]=fileparts(args.save_template_file);
                            if isempty(f1)
                                if isempty(args.output_folder)
                                    f1=pwd;
                                else
                                    f1=args.output_folder;
                                end
                                if inputDataIsCsv
                                    f3=fileparts(args.csv_file_or_data);
                                    if ~isempty(f3)
                                        f1=f3;
                                    end
                                end
                            end
                        elseif inputDataIsCsv
                            [f1, f2]=fileparts(csv_file_or_data);
                            if isempty(f1)
                                f1=pwd;
                            end
                            f2=[f2 '.umap'];
                        else
                            f1=[];
                        end
                        if ~isempty(f1) && ~isempty(pythonTemplate)
                            umap.pythonTemplate=fullfile(f1, [f2 '.python']);
                            movefile(pythonTemplate, umap.pythonTemplate, 'f');
                        end
                    end
                    if  ~isempty(args.save_template_file)
                        f1=fileparts(args.save_template_file);
                        if isempty(f1)
                            f1=pwd;
                            if inputDataIsCsv
                                f2=fileparts(args.csv_file_or_data);
                                if ~isempty(f2)
                                    f1=f2;
                                end
                            end
                            args.save_template_file...
                                =fullfile(f1, args.save_template_file);
                        end
                        args.save_template_file=...
                            File.ExpandHomeSymbol(args.save_template_file);
                        save(args.save_template_file, 'umap');
                    end
                    umap.restoreSupervisorGuiStuff(stuff, oldArgs);
                end
                if nArgOut>2 || ~strcmpi(args.cluster_output, 'none') ...
                        || ~isempty(legendRois) 
                    if nArgOut<3 ...
                            && ~strcmpi(args.cluster_output, 'graphic') ...
                            && isempty(legendRois) 
                       %if legendRois exists, UST is visble, so clustering 
                       %state must be seen too
                        warning('No clusterIdentifiers output argument');
                    elseif ~strcmpi(args.cluster_output, 'ignore')
                        try
                            clusterIdentifiers=doClusters();
                            if isempty(clusterIdentifiers)
                                dispNoDbScan;
                            end
                        catch ex
                            ex.getReport
                        end
                    end
                end
                if exist('pu', 'var') && isa(pu, 'PopUp')
                    pu.stop;
                    pu.dlg.dispose;
                end    
            end
            
            function reduceHiD2LoD
                if isSupervising
                    args.reductionType=UMAP.REDUCTION_SUPERVISED_TEMPLATE;
                else
                    if isempty(template_file)
                        if hasLabels
                            args.reductionType=UMAP.REDUCTION_SUPERVISED;
                        else
                            args.reductionType=UMAP.REDUCTION_BASIC;
                        end
                    else
                        args.reductionType=UMAP.REDUCTION_TEMPLATE;
                    end
                    if args.matchingUst
                        warning('qf_dissimilarity=true and match_scenarios=1 or 2 ONLY affects supervised template reduction')
                        args.match_scenarios(args.ustMatches)=[];
                        args.matchingUst=false;
                    end
                end
                
                if args.match_predictions 
                    if strcmp(UMAP.REDUCTION_TEMPLATE, args.reductionType)...
                            || strcmp(UMAP.REDUCTION_SUPERVISED, args.reductionType)
                        warning(['Reduction is %s ... argument match_predictions='...
                            'true applies ONLY for reductions with basic or ' ...
                            'supervised template reductions'], args.reductionType);
                    end
                    if isempty(testSetLabels)
                        warning(['Argument match_predictions=true requires label_column ' ...
                            'argument to indicate a prior classification']);
                    end
                end
                
                if args.confusion_chart
                   if strcmp(UMAP.REDUCTION_TEMPLATE, args.reductionType)...
                            || strcmp(UMAP.REDUCTION_SUPERVISED, args.reductionType)
                        warning(['Reduction is %s ... argument confusion_chart=true '...
                            'applies ONLY for reductions with basic or ' ...
                            'supervised template reductions'], args.reductionType);
                    end
                    if isempty(testSetLabels)
                        warning(['Argument confusion_chart=true requires label_column ' ...
                            'argument to indicate a prior classification']);
                    end
                end
                
                if beGraphic
                    set(fig, 'name', [get(fig, 'name') ', ' args.reductionType ' reduction']);
                end
                if strcmp(args.reductionType, UMAP.REDUCTION_SUPERVISED)
                    if ~strcmpi(args.supervised_metric, UmapUtil.CATEGORICAL)
                        if ~KnnFind.CheckDistArgs(nCols, struct(...
                                'dist_args', args.supervised_dist_args, ...
                                'metric', args.supervised_metric))
                            error('Incorrect supplementary supervised metric args (P, Cov, Scale) ');
                        end
                        umap.target_metric=args.supervised_metric;
                        umap.target_metric_kwds=args.supervised_dist_args;
                    end
                else
                    if argued.contains('supervised_metric')
                        warning(['''supervised_metric'' only affects supervised '...
                            'reduction ... need label_column argument']);
                    end
                    if argued.contains('supervised_dist_args')
                        warning(['''supervised_dist_args'''' only affects supervised '...
                            'reduction ... need label_column argument']);
                    end
                end
                extras.args=args;
                strReduction=[UmapUtil.GetReductionLongText(args.reductionType) ' reduction'];
                if args.fast_approximation
                    reportProgress(['Approximating ' strReduction ', v' UMAP.VERSION], true);
                else
                    reportProgress(['Running ' strReduction ', v' UMAP.VERSION], true);
                end
                if ~beQuiet
                    disp(UMAP.REDUCTION_TYPES);
                end
                if ~isempty(umap.supervisors)
                    umap.supervisors.setGraphicsArgs(args);
                    % semi-hack below since not sure if some developers depend on the
                    % original input version of args.parameter_names that MIGHT differ
                    % from parameter_names.  Supervisors object however needs the final
                    % version of parameter_names
                    umap.supervisors.setArgs('parameter_names', parameter_names);
                end
                if ~isempty(umap.supervisors)
                    densityBars=DensityBars([inData;...
                        umap.raw_data]);
                else
                    densityBars=DensityBars(inData);
                end
                if args.fast_approximation
                    reportProgress('Probability binning for fast approximation', true);
                    if beGraphic
                        pu.setText2('Compressing for fast approximation...');
                    end
                    probability_bins=SuhProbabilityBins(inData, true);
                    if ~isempty(umap.supervisors)
                        umap.supervisors.probability_bins=probability_bins;
                    end
                    inData=probability_bins.compress;
                    if hasLabels
                        labels=probability_bins.fit(labels);
                    end
                    if beGraphic
                        pu.setText2(...
                            sprintf('(fast approximation condenses %s rows to %s)',...
                            String.encodeK(nRows), ...
                            String.encodeK(size(probability_bins.weights,1))));
                    end
                end
                if ~isempty(template_file)
                    if args.phate
                        msgWarning(Html.WrapHr(['Reducing with UMAP' ...
                            '<br>PHATE is not supported with' ...
                            ' templates']), 0, 'north east');
                        args.reduction_algorithm='UMAP';
                        args.save_template_file='';
                    end
                    if ~isempty(umap.pythonTemplate)
                        args.python=true;
                    end
                    if ~args.python
                        superviseByMlp=false;
                        if isSupervising && args.mlp_only
                            if ~args.mlp_supervise
                                warning(['mlp_only==true thus ' ...
                                    'setting mlp_supervise to true']);
                                args.mlp_supervise=true;
                            end
                        end
                        if isSupervising && args.mlp_supervise
                            if args.fast_approximation
                                superviseByMlp=sprv.predictMlp(probability_bins, pu);
                            else
                                superviseByMlp=sprv.predictMlp(inData, pu);
                            end
                        end
                        if superviseByMlp
                            if args.mlp_only
                                if beGraphic
                                    close(fig);
                                end
                                if args.fast_approximation
                                    sprv.resolveTestDataMatching(...
                                        probability_bins, ...
                                        args.match_supervisors, pu, ...
                                        args.mlp_confidence);
                                else
                                    sprv.resolveTestDataMatching( ...
                                        inData, ...
                                        args.match_supervisors, pu, ...
                                        args.mlp_confidence);
                                end
                                if beGraphic && isempty(args.conclude)
                                    msg(Html.WrapHr([...
                                        'No reduction requested...<br>'...
                                        'predicted MLP labels are in' ...
                                        '<br><b>umap.supervisors.' ...
                                        'mlp_labels']));
                                end
                                return;
                            end
                            reduction=umap.fit_transform(inData, sprv.mlp_supervise_labels);
                        elseif ~args.joined_transform
                            reduction=umap.transform(inData);
                        else
                            reduction=umap.transform2(inData);
                        end
                    else
                        if isempty(umap.pythonTemplate) || ~exist(umap.pythonTemplate, 'file')
                            reduction=[];
                            msg(Html.WrapHr('Python template file not found'),  8,...
                                'south', 'Error...', 'error.png');
                        else
                            inFile=[tempname '.csv'];
                            reduction=UmapPython.Go(inData,inFile, [], ...
                                lower(umap.metric), umap.n_neighbors, ...
                                umap.min_dist, umap.n_components, [], ...
                                umap.pythonTemplate, args.verbose);
                        end
                    end
                else
                    if args.phate && ...
                        (~isempty(args.save_template_file) || hasLabels)
                        msgWarning(Html.WrapHr(['Reducing with UMAP' ...
                            '<br>PHATE is not supported with' ...
                            '<br>templates or supervision']), 0, ...
                            'north east');
                        args.reduction_algorithm='UMAP';
                        args.save_template_file='';
                    end
                    if ~isempty(args.mlp_train)
                        if ~isempty(args.save_template_file)
                            [pth, fn, ext]=fileparts(args.save_template_file);
                            if isempty(pth)
                                %look like a run_umap/examples case
                                if args.fast_approximation
                                    args.save_template_file=...
                                        [fn UmapUtil.FAST_FILE_SUFFIX ext];
                                end
                            end
                        end
                        if isempty(args.save_template_file) || ~hasLabels
                            warning('mlp_train argument is ONLY used when creating a supervised template');
                        else
                            [usePython, mlpArgs]=UmapUtil.GetMlpTrainArg(args.mlp_train);
                            f1=fileparts(args.save_template_file);
                            if isempty(f1)
                                f1=pwd;
                                if inputDataIsCsv
                                    f2=fileparts(args.csv_file_or_data);
                                    if ~isempty(f2)
                                        f1=f2;
                                    end
                                end
                                args.save_template_file...
                                    =fullfile(f1, args.save_template_file);
                            end
                            mlpArgs{end+1}='model_default_folder';
                            mlpArgs{end+1}=f1;
                            mf=Args.Get('model_file', mlpArgs{:});
                            if isempty(mf)
                                mlpArgs{end+1}='model_file';
                                mlpArgs{end+1}=...
                                    Supervisors.MlpFile(...
                                    args.save_template_file);
                            end
                            if exist('pu', 'var') && ~isempty(pu)
                                drawnow;
                                mlpArgs{end+1}='pu';
                                mlpArgs{end+1}=pu;
                                if ~isempty(pu.cancelBtn)
                                    pu.cancelBtn.setVisible(false);
                                end
                            end
                            try
                                mlpArgs{end+1}='label_file';
                                mlpArgs{end+1}=args.label_file;
                                umap.train_mlp([inData labels],...
                                    usePython, mlpArgs{:});
                            catch ex
                                halt=~askYesOrNo(Html.WrapHr(...
                                    ['MLP training errors...' ...
                                    '<br><br>So ...continue <i>without</i>' ...
                                    ' MLP training?']));
                                Gui.MsgException(ex);
                                if halt
                                    throw(ex);
                                end
                            end
                            if exist('pu', 'var') && ~isempty(pu)
                                pu.setText('Supervising UMAP with MLP model');
                                pu.dlg.pack;
                                %pu.dlg.setSize(puSize);
                                if ~isempty(pu.cancelBtn)
                                    pu.cancelBtn.setVisible(true);
                                end
                            end
                        end
                    end
                    if ~args.python
                        if ~hasLabels
                            if args.phate   
                                if beGraphic
                                    pu.setText2('Running PHATE from the Krishnaswamy Lab');
                                    cb=pu.addListToSouth(1.7);
                                else
                                    cb=args.verbose;
                                end
                                if isempty(args.args_phate)
                                    reduction=phate(inData, ...
                                        'k', args.n_neighbors,...
                                        'distfun', args.metric,...
                                        'n_landmarks', args.landmarks,...
                                        'ndim', args.n_components, ...
                                        'progress_callback', cb);
                                else
                                    if isempty(args.args_phate.progress_callback)
                                        args.args_phate.progress_callback...
                                            =cb;
                                    end
                                    varArg=SuhStruct.ToNamedValueCell( ...
                                        args.args_phate);
                                    reduction=phate(inData, varArg{:});
                                end
                            else
                                reduction=umap.fit_transform(inData);
                            end
                        else
                            reduction=umap.fit_transform(inData, labels);
                            if ~isempty(reduction)
                                if ~isempty(labelMap)
                                    umap.setSupervisors(labels, ...
                                        labelMap, curAxes, args);
                                end
                            end
                        end
                    else
                        inFile=[tempname '.csv'];
                        if ~hasLabels
                            labels=[];
                        end
                        reduction=UmapPython.Go(inData,inFile, [], ...
                            lower(umap.metric), umap.n_neighbors, ...
                            umap.min_dist, umap.n_components, labels, ...
                            [], args.verbose);
                        pythonTemplate=fullfile(fileparts(inFile), ...
                            [UmapPython.PYTHON_TEMPLATE '.umap']);
                        if ~isempty(reduction)
                            umap.embedding=reduction;
                            umap.raw_data=inData;
                            if ~isempty(labelMap)
                                umap.setSupervisors(labels, ...
                                    labelMap, curAxes, args);
                            end
                        end
                    end
                    if ~isempty(umap.supervisors)
                        sprv=umap.supervisors;
                        sprv.description=args.description;
                    end
                end
                if args.rescale~=1
                    reduction=reduction*args.rescale;
                end
                if args.fast_approximation && ~isempty(reduction)
                    if isSupervising
                        sprv.inputData=umap.raw_data;
                        sprv.resolveTestDataMatching(probability_bins, ...
                            args.match_supervisors, pu, ...
                            args.mlp_confidence);
                    end
                    labels=probability_bins.originalLabels;
                    inData=probability_bins.originalData;
                    reduction=probability_bins.decompress(...
                        reduction, args.jiggle_probability_bin);
                elseif isSupervising
                    sprv.inputData=umap.raw_data;
                    sprv.resolveTestDataMatching( ...
                        inData, ...
                        args.match_supervisors, ...
                        pu, ...
                        args.mlp_confidence);
                end
                if ~isempty(paramAnnotation)
                    try
                        set(paramAnnotation, 'visible', 'on');
                    catch
                    end
                end
                if ~isempty(reduction)
                    mins=min(reduction);
                    maxs=max(reduction);
                end
            end
            
            
            function ok=handleBadData
                ok=false;
                if ischar(inData)
                    msgError(Html.WrapHr(['Sorry...cannot proceed<br>'...
                        '<br>(the data is invalid... SIGH!)']));
                    globals.save;
                    return;
                end
                nanRows=any(isnan(inData'));
                badRows=sum(nanRows);
                if badRows>0
                    if beGraphic
                        if askYesOrNo(Html.WrapHr(['Data matrix has ' ...
                                String.Pluralize2('row', badRows) ...
                                'with NAN values <br>which cause odd effects on '...
                                args.reduction_algorithm '!<br>Try to remove nan values?']))
                            inData=inData(~nanRows,:);
                        end
                    else
                        warning(['Data matrix has ' ...
                            String.Pluralize2('row', badRows) ...
                            'with NAN values!']);
                        inData=inData(~nanRows,:);
                    end
                    if any(isnan(inData(:)))
                        showMsg(Html.WrapHr(['Sorry...<br>cannot proceed<br>'...
                            '<br>NAN values exist... SIGH!']));
                        globals.save;
                        return;
                    end
                end
                args.hiD=nCols-labelCols;
                if strcmpi(method, 'C++')
                    if ~StochasticGradientDescent.IsAvailable
                        if ~askYesOrNo(Html.Wrap(...
                                ['This C++ executable is missing or corrupt:'...
                                '<br>"<b>' StochasticGradientDescent.GetCmd '</b>"'...
                                '<br><br>Maybe try rebuilding by changing Clang++ '...
                                'to g++ in the build scripts in the same folder...<br>'...
                                '<br><center>Try <b>method=Java</b> instead?</center><hr>']))
                            return;
                        end
                        method='MEX';
                    end
                end                
                ok=true;
            end
            
            function ok=prepareUmapSettings
                ok=false;
                umap.init=args.init;
                isSupervising=isprop(umap, 'supervisors') && ~isempty(umap.supervisors);
                if args.matchingTestLabels
                    if isempty(template_file)
                        if length(testSetLabels) ~= length(labels)
                            testSetLabels=labels;
                        end
                    end
                end                
                if isSupervising
                    if args.n_components==2
                        progressMatchType=0;
                    else
                        limit=args.match_3D_limit;
                        if limit<nRows
                            progressMatchType=-1;
                        else
                            progressMatchType=3;
                        end
                    end
                else
                    if any(args.match_supervisors>0)
                        if argued.contains('match_supervisors')
                            if ~isequal(1, args.match_supervisors)
                                warning('match_supervisors only affects supervised template reduction');
                                args.match_supervisors=1;
                            end
                        end
                    end
                end
                UmapUtil.SetArgsTemplateCanOverride(umap, args, argued, parameter_names);
                umap.n_epochs=args.n_epochs;
                umap.nn_descent_min_rows=args.nn_descent_min_rows;
                umap.nn_descent_min_cols=args.nn_descent_min_cols;
                umap.nn_descent_max_neighbors=args.nn_descent_max_neighbors;
                umap.nn_descent_transform_queue_size=args.nn_descent_transform_queue_size;
                umap.eigen_limit=args.eigen_limit;
                umap.probability_bin_limit=args.probability_bin_limit;
                if strcmpi('Java', args.method)
                    if ~initJava
                        args.method='MEX';
                        showMsg(Html.WrapHr('Could not load umap.jar for Java method'), ...
                            'Problem with Java...', 'south west', false, false);
                    end
                elseif strcmpi('Mex', args.method)
                    exeSGD=fullfile(curPath, UmapUtil.LocateMex('sgd'));
                    exeNN=fullfile(curPath, UmapUtil.LocateMex);
                    if ~exist(exeSGD, 'file') || ~exist(exeNN, 'file')
                        UmapUtil.OfferFullDistribution(true)
                        globals.save;
                        if ~exist(exeSGD, 'file') || ~exist(exeNN, 'file')
                            if ~askYesOrNo(['<html>Continue more slowly <br>'...
                                    'without accelerants?<hr></html>'])
                                if beGraphic
                                    close(fig);
                                end
                                return;
                            end
                        end
                    end
                end
                method=umap.setMethod(args.method);
                umap.verbose=~beQuiet;
                umap.randomize=args.randomize;
                umap.setParallelTasks(args);
                if ~isSupervising ...
                        && argued.contains('see_training')
                    if ~isempty(template_file)
                        warning(['''see_training''==' String.toString(args.see_training) ...
                            ' has no effect because your template lacks supervisory labels.']);
                    else
                        warning(['''see_training''==' String.toString(args.see_training) ...
                            ' has no effect because you are not using a supervised template.']);
                    end
                end    
                ok=true;
            end
            
            function [ok,template_file]=prepareData
                ok=false;
                template_file=args.template_file;
                args=UmapUtil.RelocateExamples(args);
                csv_file_or_data=args.csv_file_or_data;
                if nArgOut>=3 && args.n_components>2
                    % check for presence of DBSCAN
                    if ~Density.HasDbScan(false)
                        if beGraphic
                            if ~askYesOrNo(Html.WrapHr(['DBSCAN for clustering in 3+D is '...
                                    '<br>not downloaded ...Continue?']))
                                delete(fig);
                                globals.save;
                                return;
                            end
                        end
                        dispNoDbScan;
                    end
                end
                if inputDataIsCsv
                    if ~exist(csv_file_or_data, 'file')
                        if beGraphic
                            delete(fig);
                        end
                        if startsWith(csv_file_or_data, UmapUtil.LocalSamplesFolder)
                            [~,fn, ext]=fileparts(csv_file_or_data);
                            if askYesOrNo(struct('icon', 'error.png', ...
                                    'msg', ['<html>Cannot access the example "' ...
                                    fn ext '"<br><br>' globals.h2Start '<center>Open '...
                                    'our shared Google Drive?</center>' globals.h2End '<hr></html>']),...
                                    'Example not found')
                                UmapUtil.GoogleDrive([], false, true);
                                msg('All sample data is in examples subfolder');
                            end
                        else
                            msg(['<html>The CSV file <br>"<b>' ...
                                globals.smallStart csv_file_or_data globals.smallEnd ...
                                '</b>"<br><center><font color="red"><i>cannot be found!!' ...
                                '</i></font><hr></center></html>'], 25, 'center', ...
                                'Error...', 'error.png');
                        end
                        globals.save;
                        return;
                    end
                    if beGraphic
                        dd=dir(csv_file_or_data);
                        set(get(fig, 'currentAxes'), 'XTick', [], 'YTick', []);
                        if dd.bytes>10*1000000 % 10 million
                            loadingHtml=['<html>' args.reduction_algorithm ' is loading '...
                                '<b><i>' String.encodeGb(dd.bytes,[],2) ...
                                '</i></b> from' Html.FileTree(csv_file_or_data) ...
                                '<br><br>'...1
                                Gui.UnderConstructionImg(.7)  '<hr></html>'];
                            an=Gui.ImageLabel(loadingHtml, 'umap.png', ...
                                '',[],fig);
                            pos=get(an,'position');
                            set(an, 'position', [100 100 pos(3) pos(4)]);
                        else
                            an=[];
                        end
                        drawnow;
                    end
                    [inData, parameter_names]=File.ReadCsv(csv_file_or_data);
                    if beGraphic
                        if ~isempty(an)
                            delete(an);
                        end
                    end
                else
                    if beGraphic
                        set(get(fig, 'currentAxes'), 'XTick', [], 'YTick', []);
                        drawnow;
                    end
                    inData=csv_file_or_data;
                    parameter_names=args.parameter_names;
                end
                % fill out args for callbacks like save_roi
                args.unreduced_data=inData; 
                if isempty(parameter_names)
                    C=size(inData,2);
                    args.unreduced_column_names=cell(1,C);
                    for c=1:C
                        args.unreduced_column_names{c}=['Column ' num2str(c)];
                    end
                else
                    args.unreduced_column_names=parameter_names;
                end
                newSubsetIdxs=[];
                template_file=args.template_file;
                [nRows, nCols]=size(inData);
                if ~KnnFind.CheckDistArgs(nCols, args)
                    error('Incorrect supplementary metric args (P, Cov, Scale) ');
                end                
                sCols=num2str(nCols);
                if nRows*nCols<15
                    if isscalar(inData)
                        if inData>1&&inData<=34
                            if askYesOrNo(['Run example #' num2str(inData) '??'])
                                run_examples(inData)
                                return
                            end
                        end
                    end
                    exitNow=true;
                    if isempty(inData)
                        if beGraphic
                            delete(fig);
                        end                        
                        if inputDataIsCsv && startsWith(csv_file_or_data, ...
                                UmapUtil.LocalSamplesFolder) && ...
                                askYesOrNo(Html.WrapHr('Remove corrupt example file?'))
                            tempfile=[tempname '.html'];
                            movefile(csv_file_or_data, tempfile);
                            if askYesOrNo(['<html>' Html.H2('File removed ...')...
                                    '<br>Inspect the contents in your browser?</html>'])
                                Html.BrowseFile(tempfile);
                            end
                        else
                            msgError('Table of numbers NOT provided');
                        end
                    else
                        if ~strcmpi(args.verbose, 'none')
                            answ=questDlg(struct('icon', 'error.png', 'msg', sprintf([...
                                '<html>Data seems too little: %d row(s) X  %d col(s)...'...
                                '<br><br><center>Proceed with this very small'...
                                ' reduction??</center><hr></html>'],...
                                nRows, nCols)), 'Wow ... not much to do?');
                            exitNow=~strcmpi('yes', answ);
                            if exitNow && beGraphic
                                delete(fig);
                            end
                        end
                    end
                    if exitNow
                        return
                    end
                end
                if strcmpi(args.label_column, 'end')
                    args.label_column=nCols;
                end
                if args.label_column>0
                    if  args.label_column>nCols
                        msg(Html.WrapHr(['The input data has ' sCols ' columns ...<br>'...
                            'THUS the label_column must be >=1 and <= ' sCols]));
                        assert(args.label_column>0 && args.label_column<=nCols, [...
                            'label_column must be >=1 and <= ' sCols]);
                    end
                    labelCols=1;
                    if args.matchingTestLabels
                        testSetLabels=inData(:, args.label_column);
                    end
                    if ~isempty(template_file)
                        if args.label_column<=length(parameter_names)
                            parameter_names(args.label_column)=[];
                        end
                    else
                        labels=inData(:, args.label_column);
                    end
                    inData(:, args.label_column)=[];
                else
                    labelCols=0;
                end
                if ~isnan(args.std_outliers)
                    inData=MatBasics.FitByStd(inData, args.std_outliers);
                end
                if strcmpi(args.metric, 'precomputed') && (~issymmetric(inData) ...
                        || ~all(diag(inData) == 0) || ~all(all(inData >=0)))
                    errTxt=['''metric''==''precomputed'', requires input data<br> '...
                        'to be a square distance matrix.'];
                    if exist('fig', 'var')
                        delete(fig);
                        msg(Html.WrapHr(errTxt), 10, 'center', ...
                            'Precomputed??', 'error.png');
                        drawnow;
                    end
                    error(errTxt);
                end

                if ~isempty(args.compress)
                    if labelCols==0
                        inData=LabelBasics.Compress(inData, [], args.compress);
                    else
                        if isempty(template_file)
                            if isscalar(args.compress)
                                [inData,labels]=LabelBasics.Compress(...
                                    inData, labels, args.compress);
                            else
                                [inData,labels]=LabelBasics.Compress(inData, labels, ...
                                    args.compress(1), args.compress(2));
                            end
                        else
                            if isscalar(args.compress)
                                [inData,testSetLabels]=LabelBasics.Compress(inData, ...
                                    testSetLabels,args.compress);
                            else
                                [inData,testSetLabels]=LabelBasics.Compress(inData, ...
                                    testSetLabels, args.compress(1), args.compress(2));
                            end
                        end
                    end
                end
                if ~isempty(args.synthesize)
                    if labelCols==0
                        inData=UmapUtil.Synthesize(inData,[],0,args.synthesize(1), ~beQuiet);
                    else
                        if isempty(template_file)
                            if length(args.synthesize)<2
                                [inData, labels]=UmapUtil.Synthesize(inData, labels,...
                                    0, args.synthesize(1), ~beQuiet);
                            else
                                [inData, labels]=UmapUtil.Synthesize(inData, labels,...
                                    args.synthesize(2), args.synthesize(1), ~beQuiet);
                            end
                        else
                            if length(args.synthesize)<2
                                [inData, testSetLabels]=UmapUtil.Synthesize(inData,...
                                    testSetLabels, 0, args.synthesize(1), ~beQuiet);
                            else
                                [inData, testSetLabels]=UmapUtil.Synthesize(inData, labels,...
                                    args.synthesize(2), args.synthesize(1), ~beQuiet);
                            end
                        end
                    end
                end
                nRows=size(inData, 1);                
                if ~isempty(template_file)
                    if ischar(template_file)
                        if ~exist(template_file, 'file')
                            if beGraphic
                                delete(fig);
                            end
                            msg(['<html>The template file <br>"<b>' ...
                                globals.smallStart template_file globals.smallEnd ...
                                '</b>"<br><center><font color="red"><i>cannot be found!!' ...
                                '</i></font><hr></center></html>'], 25, 'center', ...
                                'Error...', 'error.png');
                            globals.save;
                            return;
                        end
                        if ~isempty(parameter_names) && length(parameter_names)~=size(inData, 2)
                            msgBox( struct( 'icon', 'error.png', 'msg', ...
                                Html.WrapHr(sprintf(['<b>Cannot create '...
                                'or use template</b> ...<br>'...
                                '%d parameter_names... but data has %d parameters?'], ...
                                length(parameter_names), size(inData,2)))));
                            if beGraphic
                                delete(fig)
                            end
                            globals.save;
                            return;
                        end
                        [umap, ~, canLoad, reOrgData, paramIdxs]=...
                            Template.Get(inData, parameter_names, ...
                            template_file, 3);
                        if ~isempty(umap)
                            umap.args=args;
                        end
                        if ~isempty(reOrgData)
                            % column label order differed
                            inData=reOrgData;
                            nCols=size(inData, 2);
                            sCols=num2str(nCols);
                            if ~isempty(parameter_names) && ~isempty(paramIdxs)
                                parameter_names=parameter_names(paramIdxs);
                            end
                        end
                    elseif isa(template_file, 'UMAP')
                        umap=template_file;
                        if ~isempty(umap)
                            umap.args=args;
                        end
                        canLoad = true;
                    end
                    if isempty(umap)
                        if ~canLoad
                            if beGraphic
                                showMsg(Html.WrapHr(['No template data found in <br>"<b>', ...
                                    template_file '</b>"']));
                            else
                                disp(['No template data found in ' template_file]);
                            end
                        end
                        if beGraphic
                            delete(fig);
                        end
                        globals.save;
                        return;
                    else
                        args.n_components=umap.n_components;
                        if ~isempty(umap.supervisors)
                            sprv=umap.supervisors;
                            sprv.verbose=args.verbose;
                            if args.color_defaults
                                sprv.overrideColors(args.color_file, beQuiet);
                            end
                            if beQuiet
                                LabelBasics.Frequency(sprv.labels, sprv.labelMap, true)
                            else
                                LabelBasics.Frequency(sprv.labels, sprv.labelMap, true, [])
                            end
                            
                            sprv.description=args.description;
                            sprv.context=args.context;
                            
                            %Connor's NEW joined_transform immunizes reduction from
                            %false positives if items in the test set are TOO different
                            %from the training set
                            
                            if ~args.joined_transform
                                [percNewSubsets, unknownIdxs]=...
                                    Template.CheckForUntrainedFalsePositives(umap, inData);
                                if percNewSubsets>13 && beGraphic
                                    [choice, cancelled]=Template.Ask(percNewSubsets);
                                    if cancelled
                                        if beGraphic
                                            delete(fig);
                                        end
                                        globals.save;
                                        return;
                                    end
                                    if choice==2
                                        umap.clearLimits;
                                        newSubsetIdxs=unknownIdxs;
                                        template_file=[];
                                    end
                                end
                            end
                            sprv.initClustering(args.cluster_detail{1}, ...
                                args.cluster_method_2D, args.minpts, ...
                                args.epsilon, args.dbscan_distance);
                            sprv.initPlots(args.contour_percent);
                        end
                    end
                    if argued.contains('probability_bin_limit')
                        warning('Your probability_bin_limit setting has no effect when templates are used...');
                    end
                    if argued.contains('eigen_limit')
                        warning('Your eigen_limit setting has no effect when templates are used...');
                    end
                    if argued.contains('init')
                        warning('Your init setting has no effect when templates are used...');
                    end
                else
                    umap=UMAP;
                    umap.args=args;
                end
                ok=true;
            end
            
            function prepareLabels
                nParams=length(parameter_names);
                good=nParams==0||nParams==nCols || (args.label_column>0 &&...
                    (nParams==nCols-1 || nParams==nCols));
                if ~good
                    if args.label_column>0
                        preAmble=sprintf(['# data columns=%d, # parameter_names=%d '...
                            'since label_column=%d <br># parameter_names must be '...
                            '%d or %d '],  nCols, nParams, args.label_column, ...
                            nCols, nCols-1);
                    else
                        preAmble=sprintf(['# of data columns(%s) must equal '...
                            '# of parameter_names(%d)'], sCols, nParams);
                    end
                    msg(Html.WrapHr(preAmble));
                    assert(nParams==0||nParams==nCols || (args.label_column>0 &&...
                        (nParams==nCols-1 || nParams==nCols)), preAmble);
                end
                if ~isempty(newSubsetIdxs)
                    hasLabels=true;
                    labelCols=0;
                    [labels, labelMap]=resupervise(umap, inData, newSubsetIdxs);
                    nLabels=length(unique(labels));
                elseif args.label_column>0 && isempty(template_file) && ~args.matchingUmap
                    hasLabels=true;
                    labelCols=1;
                    good=args.label_column>0 && args.label_column<=nCols;
                    if ~good
                        msg(Html.WrapHr(['The input data has ' sCols ' columns ...<br>'...
                            'THUS the label_column must be >=1 and <= ' sCols]));
                        assert(args.label_column>0 && args.label_column<=nCols, [...
                            'label_column must be >=1 and <= ' sCols]);
                    end
                    nLabels=length(unique(labels));
                    if nLabels > .2*nRows
                        if ~acceptTooManyLabels(nLabels)
                            if beGraphic
                                delete(fig);
                            end
                            return;
                        end
                        preAmble='%d is a LOT of distinct labels for a %dx%d matrix!';
                        msg(['WARNING:  ' sprintf(preAmble, nLabels, nRows, nCols)]);
                        warning(preAmble, nLabels, nRows, nCols);
                    end
                    if args.label_column<=nParams
                        parameter_names(args.label_column)=[];
                    end
                    umap.dimNames=parameter_names;
                    nLabels=length(unique(labels));
                else
                    hasLabels=false;
                    nLabels=0;
                end
                if args.label_column>0
                    if isSupervising && ~isempty(args.label_file) && isempty(testSetLabels)
                        if args.label_column==0
                            warning(['label_file has no effect when reducing '...
                                'with supervised template without test set labels']);
                        elseif ~any(args.match_scenarios==1) ....
                                && ~any(args.match_scenarios>2)
                            warning(['test set labels only needed for UMAP'...
                                ' supervised templates IF specifying '...
                                'match_scenarios 1, 3, or 4 ']);
                        end
                    else
                        if ~isempty(testSetLabels)
                            [labelMap, halt]=getLabelMap(testSetLabels);
                        else
                            [labelMap, halt]=getLabelMap(labels);
                            if hasLabels
                                ColorsByName.Override(labelMap, args.color_file, beQuiet);
                            end
                        end
                        if halt
                            if beGraphic
                                delete(fig);
                            end
                            globals.save;
                            return;
                        end
                        if ~isempty(labelMap) && ~args.buildLabelMap ...
                                && isSupervising && ~isempty(testSetLabels)
                            testSetLabels=LabelBasics.RelabelIfNeedBe( ...
                                testSetLabels, sprv.labelMap, labelMap);
                        end
                        
                    end
                end
            end
            
            function prepareAnnotations
                info=[String.encodeK(nRows) ' x ' String.encodeInteger(nCols-labelCols)];
                if inputDataIsCsv
                    [~, fileName]=fileparts(csv_file_or_data);
                    figName=[args.reduction_algorithm ' on ' fileName '. ' info];
                    info=[args.reduction_algorithm ' on ' String.ToLaTex(fileName) ', ' info];
                else
                    info=[args.reduction_algorithm ' on ' info];
                    figName=info;
                end
                if args.python
                    info=[info ', Python'];
                else
                    info=[info ', ' method];
                end
                set(fig, 'NumberTitle', 'off', 'Name', figName );
                drawnow;                
                pause(.01);
                if strcmpi(method, 'Java') || strcmpi(method, 'C++') || strcmpi(method, 'MEX')
                    pu=initProgressPopUp;
                end
                tic;                
                info2=['(optimize\_layout method=' method ')'];
                busy=Gui.ShowBusy(fig, '<br><br><br><br><br><br>', ...
                    'simplicialComplex.png', simplicialSize);
                if ispc
                    runAnnotation=Gui.TextBox({['\color{blue}Running '...
                        info], ['\fontsize{9}' info2]}, fig, 'visible', 'off');
                else
                    runAnnotation=Gui.TextBox({['\color{blue}Running '...
                        info], ['\fontsize{11}' info2]}, fig, 'visible', 'off');
                end
                lblP=get(runAnnotation, 'position');
                lblP(2)=.008;
                set(runAnnotation, 'visible', 'on', 'position', lblP);
                updatePlot;
                strMetric=umap.metric;
                if ~ischar(strMetric)
                    strMetric='custom';
                end
                if ~isempty(umap.dist_args)
                    if strcmpi(strMetric,'minkowski')
                        strMetric=[strMetric ' P=' String.encodeBank(umap.dist_args)];
                    elseif strcmpi(strMetric,'mahalanobis')
                        strMetric=[strMetric ' Cov'];
                    elseif strcmpi(strMetric, 'seuclidean')
                        strMetric=[strMetric ' Scale'];
                    end
                end
                if ~isempty(args.args_phate)
                    txt=sprintf(['\\color{black}' ...
                        'k=\\color{blue}%d\\color{black}, '...
                        'n\\_landmarks=\\color{blue}%s\\color{black}, '...
                        'distfun=\\color{blue}%s\\color{black}, '...
                        'a=\\color{blue}%s\\color{black}, '...
                        'labels=\\color{blue}%d'], ...
                        args.args_phate.k, ...
                        num2str(args.args_phate.n_landmarks), ...
                        args.args_phate.distfun,...
                        num2str(args.args_phate.a),...
                        nLabels);
                elseif umap.spread == 1
                    txt=sprintf(['\\color{black}n\\_neighbors=\\color{blue}%d\\color{black}, '...
                        'min\\_dist=\\color{blue}%s\\color{black}, '...
                        'metric=\\color{blue}%s\\color{black}, '...
                        'randomize=\\color{blue}%d\\color{black}, '...
                        'labels=\\color{blue}%d'], ...
                        umap.n_neighbors, num2str(umap.min_dist), strMetric,...
                        umap.randomize, nLabels);
                else
                    txt=sprintf(['\\color{black}n\\_neighbors=\\color{blue}%d\\color{black}, '...
                        'min\\_dist=\\color{blue}%s\\color{black}, '...
                        'spread=\\color{blue}%s\\color{black}, '...
                        'metric=\\color{blue}%s\\color{black}, '...
                        'randomize=\\color{blue}%d\\color{black}, '...
                        'labels=\\color{blue}%d'], ...
                        umap.n_neighbors, num2str(umap.min_dist), num2str(umap.spread),...
                        strMetric, umap.randomize, nLabels);
                end
                paramAnnotation=Gui.TextBox(txt, fig, 'visible', 'off');
                pp=get(paramAnnotation, 'position');
                pp(2)=.94;
                set(paramAnnotation, 'position', pp, 'visible', 'on');
                drawnow;
            end
            
            function progressPopUp=initProgressPopUp
                if isempty(args.progress_callback)
                    umap.progress_callback=@(javaObject)progress_report(javaObject);
                else
                    umap.progress_callback=args.progress_callback;
                end
                try
                    nTh=edu.stanford.facs.swing.StochasticGradientDescent.EPOCH_REPORTS+3;
                    figure(fig);
                    if args.qf_tree
                        puLocation='north++';
                    else
                        if args.see_training
                            puLocation='south west+';
                        else
                            puLocation='south+';
                        end
                    end
                    if args.fast_approximation ...
                            && nRows*nCols<UmapUtil.MINIMUM_FAST_APPROXIMATION
                        warning(['Fast approximation ignored for %s data points'...
                            '( minimum is %s )'], num2str(nRows*nCols), ...
                            num2str(UmapUtil.MINIMUM_FAST_APPROXIMATION));
                        args.fast_approximation=false;
                    end
                
                    if args.fast_approximation
                        ttl='Fast parameter reduction....';
                    else
                        ttl='Parameter reduction...';
                    end
                    imgFile=Gui.GetResizedImageFile( ...
                        'simplicialComplex.png', .1);
                    progressPopUp=PopUp(Html.WrapHr(sprintf([args.reduction_algorithm ...
                        ' is reducing  <b>%d</b> parameters to ' ...
                        num2str(args.n_components) '...'], nCols-labelCols)), ...
                        puLocation, ttl, false, true, imgFile, false, ...
                        [], Gui.JWindow(fig), true);
                    progressPopUp.initProgress(nTh);
                    progressPopUp.pb.setStringPainted(true);
                    progressPopUp.setTimeSpentTic;
                    SuhWindow.Follow(progressPopUp.dlg, fig, puLocation, true);
                    progressPopUp.dlg.setVisible(true);
                    drawnow;
                catch
                    args.method='MEX';
                    method=umap.setMethod(args.method);
                    showMsg(Html.WrapHr(['Could not load umap.jar for Java method'...
                        '<br><br>Switching optimize_layout method to "MEX" ']), ...
                        'Problem with Java...', 'south west', false, false);
                end
            end
            
            function hideAnnotations
                try
                    delete(paramAnnotation);
                    if args.hide_reduction_time
                        delete(timeAnnotation);
                    end
                catch
                end
            end
            
            function [name, named]=nameTheRegion(theRoi, ask)
                name='';
                named=false;
                if nargin<2
                    ask=true;
                    if nargin<1
                        theRoi=roi;
                    end
                end
                if isempty(theRoi)
                    msg(Html.WrapHr(...
                    'No current region<br>of interest!'), ...
                    7, 'south east+');
                else
                    [txt, names, ~]=getTrainers(theRoi);
                    if ~isempty(names) && ~any(contains(names, 'no subset'))
                        if ~isempty(txt)
                            tip=Html.Wrap(txt, 200);
                            shakeNameBtn(true, tip);
                        else
                            tip=[];
                        end
                        named=true;
                        if ask
                            name=LabelLegendRoi.Rename(roiMap, theRoi, names{1}, tip);
                        else
                            name=names{1};
                        end
                    else
                        mdn=median(reduction(RoiUtil.GetRows( ...
                            theRoi, reduction), :));
                        name=sprintf('%s X=%s Y=%s', ...
                            args.reduction_algorithm,...
                            String.encodeRounded(mdn(1),1), ...
                            String.encodeRounded(mdn(2), 1));
                        if ask
                            name=LabelLegendRoi.Rename(roiMap, theRoi, name);
                        end
                    end
                end
            end
            
            
            function saveNamedRegions(name)
                dfltFile=['umap_' args.reductionType...
                    '_roi' args.output_suffix '.properties' ];
                if LabelLegendRoi.Save(roiMap, args.save_output, ...
                        args.output_folder, dfltFile, ...
                        askedToSeeOutputFolder, reduction, ...
                        args, name, sprv)
                    askedToSeeOutputFolder=true;
                elseif isempty(args.save_roi)
                    shakeNameBtn;
                end
            end
            
            function shakeNameBtn(showTip, tip)
                if isempty(btnNameRoi)
                    return
                end
                edu.stanford.facs.swing.Basics.Shake(btnNameRoi,5);
                if nargin>0 && showTip
                    if nargin<2
                        tip=Html.WrapHr([...
                        'This region of interest does not yet have<br>'...
                        'a name ... click here to give it one!']);
                    end
                    globals.showToolTip(btnNameRoi, tip, 20, 15);
                end
            end
            
            function shakeSyncBtn(showTip)
                if nargin>0 && showTip
                    globals.showToolTip(btnNameRoi, Html.WrapHr([...
                        'Click <b>here</b> to keep <br>'...
                        'sychronizing the ParameterExplorer!']), 20, 15);
                    edu.stanford.facs.swing.Basics.Shake(cbSyncKld, 8);
                else
                    edu.stanford.facs.swing.Basics.Shake(cbSyncKld, 5);
                end
            end
            function setClusterDetail
                if isempty(comboClu) || ~exist('dbm', 'var')
                    return;
                end
                if isempty(dbm) || ~strcmpi('dbm', dbm.method)
                    startIdx=0;
                else
                    dtls=Density.DETAILS(1:end-2);
                    startIdx=StringArray.IndexOfIgnoreCase(dtls, dbm.detail);
                end
                comboClu.setSelectedIndex(startIdx);
                setClusterCount;
            end
            
            function setClusterCount
                if isempty(comboClu)
                    return;
                end
                if isempty(dbm) || ~strcmpi('dbm', dbm.method)
                    html='Cluster detail ';
                else
                    html=sprintf('%d clusters ', dbm.numClusters);
                end
                lblClu.setText(Html.WrapSm( html, globals));
            end
            
            function ok=handleNoData
                ok=false;
                set(curAxes, 'xtick', []);set(curAxes, 'yTick',[]);
                if isempty(csv_file_or_data)
                    exeSGD=fullfile(curPath, UmapUtil.LocateMex('sgd'));
                    exeNN=fullfile(curPath, UmapUtil.LocateMex);
                    if ~exist(exeSGD, 'file') || ~exist(exeNN, 'file')
                        UmapUtil.OfferFullDistribution(true)
                        globals.save;
                        if ~exist(exeSGD, 'file') || ~exist(exeNN, 'file')
                            return;
                        end
                        csv_file_or_data='sample10k.csv';
                        answer=1;
                    else
                        [answer, cancelled]=Gui.Ask(Html.Wrap(...
                            ['Options for getting examples & accelerants<br>'...
                            'from the Herzenberg Lab@Stanford University.<hr>']),...
                            {'<html>Download from Google Drive <i>and exit</i></html>', ...
                            '<html>Download examples+accelerants then <i>keep running</i></html>!!'...
                            '<html>Download accelerants <b>only</b> then <i>keep running</i></html>'}, ...
                            'runUmapDownload', 'Getting our examples...', 2 );
                        if ~cancelled
                            if answer==3
                                if UmapUtil.DownloadAdditions
                                    csv_file_or_data='sample10k.csv';
                                end
                            elseif answer==2
                                if UmapUtil.DownloadAdditions(false)
                                    csv_file_or_data=downloadCsv;
                                end
                            elseif answer==1
                                UmapUtil.GoogleDrive([], false);
                            end
                        end
                    end
                    if isempty(csv_file_or_data)
                        if beGraphic
                            delete(fig);
                        end
                        globals.save;
                        return;
                    end
                    if answer==2 && ~askYesOrNo(Html.Wrap([...
                            'Test CSV files have been downloaded:<ol>'...
                            ' <li>sample10k<li>sample30k<li>sampleBalbcLabeled55k'...
                            '<li>sample130k<li>sampleRag148k.csv<li>sampleRag55k.csv'...
                            '</ol><br><center>Run UMAP on <b>sample10k</b> now?<hr></center>']))
                        if beGraphic
                            delete(fig);
                        end
                        globals.save;
                        return;
                    end
                    args.csv_file_or_data=csv_file_or_data;
                end
                ok=true;
            end
            
            function [fig, ax, tb, btnClusterRoi, cbSyncKld, comboClu,...
                    lblClu, btnNameRoi, btnUst]=prepareGui
                if args.qf_tree
                    whereIf1stTime='south';
                elseif args.matchingUmap || args.matchingUst
                    whereIf1stTime='center';
                else
                    whereIf1stTime='onscreen';
                end
                needZoom=MatBasics.usingOlderThan('R2019a');
                if isempty(args.locate_fig)
                    if args.see_training
                        [fig, ~, personalized, location]...
                            =Gui.Figure(~needZoom, 'run_umap_fig2',[],whereIf1stTime);
                        if ~personalized
                            location(3)=location(3)*.9;
                            location(4)=location(4)*1.55;
                        end
                    else
                        [fig, ~, ~, location]=Gui.Figure(~needZoom, 'run_umap_fig',[],whereIf1stTime);
                    end
                else
                    [fig, ~, ~, location]=Gui.Figure(false);
                end
                set(fig, 'name', ['Running ' args.reduction_algorithm ...
                    '...'], 'color', 'white');
                Gui.SetToRightOf(ToolBarMethods.addSpeakerButton, fig);
                drawnow;
                cbSyncKld=[]; btnNameRoi=[];
                if isdeployed ...% fewer buttons on 1st toolbar when in MATLAB runtime
                        || ~needZoom
                    tb=ToolBar.Get(fig);
                else
                    tb=ToolBar.New(fig, false);
                end
                
                if args.n_components==2 || args.roi_table==1 || args.roi_table==3
                    if ~isempty(location)
                        set(fig, 'outerposition', location); %dont keep expanding height
                    end
                    tb.setEnabled(false);
                end
                btnUst=ToolBarMethods.addButton(tb,'ust.png', ...
                    'Manage classification via supervisors', ...
                    @(h,e)showUstOptions(h));
                btnUst.setVisible(false);
                if args.n_components==2
                    clusterCreateCallback=@(h,e)createRoiWithClusterPick();
                else
                    clusterCreateCallback=[];
                end
                if args.roi_table==1 || args.roi_table==3
                    if isempty(args.save_roi)
                        btnNameRoi=ToolBarMethods.addButton(tb, [], ...
                            'Name/save current region of interest', @(h,e)roiMenu(h),...
                            'ROI:');
                    end
                    tip=Html.WrapHr(['Use this region of interest tool<br>'...
                        'to explore each dimension''s data distribution<br>'...
                        'and Kullback-Leibler Divergence']);
                    ToolBarMethods.addButton(tb, 'ellipseGate.png', tip, ...
                        @(h,e)createRoi(h, RoiUtil.ELLIPSE));
                    ToolBarMethods.addButton(tb, 'rectangleGate.png', tip, ...
                        @(h,e)createRoi(h, RoiUtil.RECTANGLE));
                    ToolBarMethods.addButton(tb, 'polygonGate.png', ...
                        tip, @(h,e)createRoi(h, RoiUtil.POLYGON));
                    tb.jToolbar.addSeparator;
                    [btnClusterRoi, comboClu, lblClu]=Density.ComboDetail(...
                        @(h,e)clusterCallback(h), tb, ...
                        clusterCreateCallback, globals);                    
                    img=Html.ImgXy('pseudoBarHi.png', [], .819);
                    cbSyncKld=Gui.CheckBox(...
                        Html.WrapSm(['Sync ' img], globals), ...
                        globals.is('run_umap.SyncKld', true), ...
                        [], '', @(h,e)syncKld(), ...
                        ['<html>Select to synchronize region of interest<br>'...
                        'tools with the ' img ' ParameterExplorer</html>']);
                    ToolBarMethods.addComponent(tb, ...
                        Gui.FlowLeftPanelBorder(cbSyncKld));
                    tb.jToolbar.addSeparator;
                else
                    [btnClusterRoi, comboClu, lblClu]=Density.ComboDetail(...
                        @(h,e)clusterCallback(h), tb, ...
                        clusterCreateCallback, globals);
                end
                n2Tip='See HeatMap for UMAP data islands';
                Gui.AddSvgToToolBar(fig, tb)
                ToolBarMethods.addButton(tb,'heatMapHot.png', n2Tip, ...
                    @(h,e)showHeatMap)
                if isempty(args.template_file)
                    ToolBarMethods.addButton(tb, ...
                        'save16.gif', ...
                        'Save this reduction as a template', ...
                        @(h,e)saveTemplate(), Html.WrapSm('Save'));
                end

                extras.fig=fig;
                if ~isempty(args.locate_fig)
                    SuhWindow.Follow(fig, args.locate_fig);
                    SuhWindow.SetFigVisible(fig);
                else
                    Gui.FitFigToScreen(fig);
                    set(fig, 'visible', 'on');
                end                
                ax=axes('Parent', fig);
                tb.setEnabled(false)
            end
            
            function saveTemplate
                if askYesOrNo(Html.WrapHr(['Save this UMAP reduction as a<br>', ...
                        'template to accelerate reductions<br>', ...
                        'and to provide consistent view<br>'...
                        'for other compatible datasets?']), ...
                        'UMAP Templates', 'north east')

                    if length(parameter_names)~=size(inData, 2)
                        msgError(Html.WrapHr(sprintf(['<b>Cannot create '...
                            'template</b> ...<br>'...
                            '%d parameter_names ...but data has %d parameters?'], ...
                            length(parameter_names), size(inData,2))));
                    else
                        [stuff, args]=umap.prepareForTemplate(curAxes, true);
                        if inputDataIsCsv
                            Template.Save(umap, csv_file_or_data);
                        elseif ~isempty(args.template_folder)
                            Template.Save(umap, fullfile(...
                                args.template_folder, 'template.csv'));
                        elseif ~isempty(args.output_folder)
                            Template.Save(umap, fullfile(...
                                args.output_folder, 'template.csv'));
                        else
                            Template.Save(umap, fullfile(pwd, 'template.csv'));
                        end
                        umap.restoreSupervisorGuiStuff(stuff, args);
                    end
                end
            end
            function lbls=getSprvLabels(rows)
                if isempty(sprvLabels)
                    sprvLabels=sprv.getLabels(probability_bins);
                end
                lbls=sprvLabels(rows);
            end
            
            function [txt, names, cnts]=getTrainers(pRoi)
                T=findQfMatch(dbm.detail);
                try                    
                    if isSupervising
                        if nargin<1 || isempty(pRoi)
                            rows=true(1, size(sprv.embedding,1));
                        else
                            rows=RoiUtil.GetRows(pRoi, sprv.embedding);
                        end
                        if nargout>1
                            [txt, names, cnts]=...
                                MatBasics.HistCountsText(...
                                sprv.labels(rows), sprv.labelMap);
                        else
                            txt=MatBasics.HistCountsText(...
                                sprv.labels(rows), sprv.labelMap);
                        end
                    else
                        if nargin<1 || isempty(pRoi)
                            rows=true(1, size(reduction,1));
                        else
                            rows=RoiUtil.GetRows(pRoi, reduction);
                        end
                        if ~isempty(T)
                            sIdPerRow=T.qf.matchTestEvents;
                        end
                        if nargout>1
                            if ~isempty(T)
                                [txt, names, cnts]=...
                                    MatBasics.HistCountsText( ...
                                    sIdPerRow(rows, 1),...
                                    T.qf.tIds, T.qf.tNames);
                            elseif ~isempty(sprv)
                                [txt, names, cnts]=...
                                    MatBasics.HistCountsText(...
                                    getSprvLabels(rows), sprv.labelMap);
                            end
                        else
                            if ~isempty(T)
                                txt=MatBasics.HistCountsText(sIdPerRow(rows, 1),...
                                    T.qf.tIds, T.qf.tNames);
                            else
                                txt=MatBasics.HistCountsText(...
                                    getSprvLabels(rows), sprv.labelMap);
                            end
                        end
                    end
                    txt=['<b>This region matches</b><hr>' txt];
                catch
                    txt='';
                    names={};
                    cnts=[];
                end
            end
            
            function roiMenu(hObj)
                jMenu=PopUp.Menu;
                Gui.NewMenuItem(jMenu, 'Name current region of interest', ...
                    @(h,e)nameTheRegion(), 'tool_data_cursor.gif');
                Gui.NewMenuItem(jMenu, 'Save regions of interest',...
                    @(h,e)saveNamedRegions(), 'save16.gif');
                jMenu.show(hObj, 5, 25);
                jMenu.repaint;
                drawnow;
            end
            
            function [mdns, nClues, nCols]=getClusterMedians(onUnreduced)
                mdns=[];
                if isempty(clusterIdentifiers)
                    return;
                end
                clues=unique(clusterIdentifiers);
                clues=clues(clues~=0);
                nClues=length(clues);
                if onUnreduced
                    D=inData;
                else
                    D=reduction;
                end
                nCols=size(D, 2);
                mdns=zeros(nClues,nCols);
                for i=1:nClues
                    clue=clues(i);
                    d=D(clusterIdentifiers==clue, :);
                    mdns(i, :)=median(d);
                end
            end
            
            function hideMeasurement
                if ~isempty(hMeasurement)
                    try
                        delete(hMeasurement);
                    catch
                    end
                    hMeasurement=[];
                    hMeasurementIdx=0;
                end
                if ~isempty(hMeasurementLabel)
                    try
                        delete(hMeasurementLabel);
                    catch
                    end
                end
                if args.n_components>2
                    updatePlot(reduction, false, true, false);
                end
            end
            
            function showMeasurement(idx)
                if idx==hMeasurementIdx
                    hideMeasurement
                    updatePlot(reduction, true, true, false);
                    showClusterLabels(4);
                    return;
                end
                hideMeasurement;
                cs=SuhHeatMap.DefaultColorScheme;
                name=parameter_names{idx};
                idx2=String.IndexOf(name, ':');
                if idx2>0
                    name=name(1:idx2-1);
                end
                html=['<html><b>' globals.smallStart name '</b> low to high: ' ...
                    Html.ImgXy([cs '.png'], [], .33) ...
                    globals.smallEnd '</html>'];
                hMeasurementLabel=Gui.PutJavaInFig(...
                    javaObject('javax.swing.JLabel', html), fig, 7, 4);
                cla(curAxes);
                vector=inData(:,idx);
                hMeasurement=scatterByColorScheme([],vector, ...
                    curAxes, reduction);
                hMeasurementIdx=idx;
                if args.n_components==2
                    hContour=updatePlot(reduction, false, true, ...
                        false, 5, false);
                    hContour.LineColor=[.8 .8 .8];
                    hContour.LineWidth=1;
                end
                showClusterLabels(2);
            end
            
            function H=drawClusterBorder3D(clue)
                l=clusterIdentifiers==clue;
                D=reduction(l,:);
                clr2='red';
                H=plot3(curAxes, D(:,1), D(:,2), D(:,3), '.', ...
                    'markerSize', 5, 'lineStyle', 'none', ...
                        'markerEdgeColor', clr2, ...
                        'markerFaceColor', clr2);
            end
            
            function hideClusterLabels
                if ~isempty(hClusterLabels)
                    try
                        for i=1:length(hClusterLabels)
                            delete(hClusterLabels{i});
                        end
                    catch
                    end
                    hClusterLabels={};
                end
            end
            
            function showClusterLabels(secs, suffix)
                if ~beGraphic || isempty(clusterIdentifiers)
                    return;
                end
                if nargin<2
                    suffix='';
                end
                hideClusterLabels;
                cntClusterLabels=cntClusterLabels+1;
                [mdns, nClues, nCols]=getClusterMedians(false);
                for i=1:nClues
                    pos=mdns(i,:);
                    if ~isempty(suffix)
                        lbl=['' String.toLetter(i) ': ' suffix];
                    else
                        lbl=String.toLetter(i);
                    end
                    if args.n_components>2
                        fs=11;
                        mrg=4;
                    else
                        fs=9;
                        mrg=2;
                    end
                    hClusterLabels{i}=text(curAxes, ...
                        'EdgeColor', 'black', 'Margin', mrg, ...
                        'FontWeight', 'bold', 'FontSize', fs, ...
                        'Position', pos, 'String', lbl, ...
                        'BackgroundColor', [.97 .97 .90]);
                end
                if nargin>0 && secs>0
                    MatBasics.RunLater(...
                        @(h,e)hide(cntClusterLabels), secs);
                end
                
                function hide(cnt)
                    if cnt==cntClusterLabels
                        hideClusterLabels();
                    end
                end
            end
            
            function showHeatMap(reClustered)
                if nargin>0 && reClustered
                    if ~isempty(jdHeatMap)
                        if ~jdHeatMap.isVisible
                            jdHeatMap.dispose;
                            jdHeatMap=[];
                            return;
                        end
                        jdHeatMap.dispose;
                        hideMeasurement;
                    else
                        return;
                    end
                elseif ~isempty(jdHeatMap)&& jdHeatMap.isVisible
                    jdHeatMap.requestFocus;
                    return;
                end

                jdHeatMap=[];
                if isempty(dbm)
                    msg(Html.WrapHr(['Select "<b>Cluster detail</b>" '...
                        'from the <br>right side of the toolbar...']), ...
                        5, 'east');
                    edu.stanford.facs.swing.Basics.Shake(comboClu, 5);
                    return;
                end
                [cnts, ~, clues]=MatBasics.HistCounts(clusterIdentifiers);
                if isempty(cnts)
                    return;
                end
                tb.setEnabled(false);
                Gui.ShowFacs(fig, 'Building HeatMap');

                nClues=length(cnts);
                [R, C]=size(inData);
                mdns=zeros(nClues, C);
                names=cell(1,nClues);
                for i=1:nClues
                    clue=clues(i);
                    l=clusterIdentifiers==clue;
                    mdns(i,:)=median(inData(l,:));
                    if clue==0
                        names{i}='Unclustered';
                    else
                        names{i}=['Cluster ' String.toLetter(clue)];
                    end
                end
                mNames=cell(1,C);
                for c=1:C
                    pn=parameter_names{c};
                    idx=String.IndexOf(pn,':');
                    if idx>0
                        mNames{c}=pn(1:idx-1);
                    else
                        mNames{c}=pn;
                    end
                end
                rowAdvice='click for ParameterExplorer';
                columnAdvice=['Click to see the level of %s' ...
                    ' expression in each data island!'];
                [jdHeatMap, tpHeat]=SuhHeatMap.New(...
                    'measurements', mdns, 'rawMeasurements', mdns,...
                    'measurementNames', mNames, ...
                    'names', names, 'freqs', cnts/R, ...
                    'columnClickedCallback', @columnClicked,...
                    'rowClickedCallback', @rowClicked,...
                    'closeFnc', @(saved, idxs, checkBoxes)close,...
                    'ignoreScatter', args.ignoreScatter,...
                    'parentFig', fig, 'ax', curAxes,...
                    'windowTitleSuffix', [' for ' args.reduction_algorithm ...
                    ' Data Islands'],...
                    'rowClickAdvice', rowAdvice, ...
                    'cellClickAdvice', rowAdvice, ...
                    'columnClickAdvice', columnAdvice,...
                    'cellClickedCallback', @cellClicked);
                updatePlot(reduction, false, true, false);
                showClusterLabels(6);
                dlg=handle(jdHeatMap, 'CallbackProperties');
                set(dlg, 'WindowClosingCallback', @(h,e)close());
                SuhWindow.Follow(jdHeatMap, fig, 'west++')
                globals.showToolTip(tpHeat, ['<html>Click a column '...
                    'heading to see <i><u>that column''s</u></i><br>level'...
                    ' of expression in each data island!<hr></html>'], ...
                    15, 5, 11, [], false, 1);
                Gui.HideBusy(fig);
                tb.setEnabled(true);

                function close
                    jdHeatMap.dispose;
                    jdHeatMap=[];
                    hideMeasurement;
                    hideClusterLabels;
                    updatePlot(reduction, false);
                    rePlotLabels;
                end
                
                function rowClicked(~, row)
                    is2D=args.n_components==2;
                    hideMeasurement;
                    clue=clues(row);
                    if clue>0
                        showClusterLabels(2);
                        if is2D
                            if clue>length(hClusterBorders) ...
                                    || ~ishandle(hClusterBorders(clue))
                                updatePlot(reduction, false, true, false);
                                hClusterBorders=dbm.drawBorders(curAxes, ...
                                    SuhRunUmap.CLUSTER_BORDER_CLR);
                            end
                            H=hClusterBorders(clue);
                            clr=get(H, 'Color');
                            lw=get(H, 'LineWidth');
                            set(H, 'LineWidth', 4, 'Color', 'red')
                            Gui.FlashN(H, 2, .25, false);
                        else
                            Gui.FlashN(hClusterLabels{row}, 3, .25, false);
                            H=drawClusterBorder3D(clue);
                        end
                    end
                    name=names{row};
                    showRoiTable(clusterIdentifiers==clue, name);
                    if clue>0
                        Gui.FlashN(H, 4, .25, ~is2D);
                        if is2D
                            MatBasics.RunLater(@(h,e)resetBorder(...
                                H, clr, lw), 4);
                        end
                    end
                end
                
                function resetBorder(H, clr, lw)
                    try
                        set(H, 'LineWidth', lw, 'Color', clr);
                    catch
                    end
                end
                
                function cellClicked(ax, row, ~)
                    rowClicked(ax, row);
                end
                
                function columnClicked(~, column)
                    showMeasurement(column);
                    if args.n_components==2
                        hClusterBorders=dbm.drawBorders(curAxes, ...
                                    SuhRunUmap.CLUSTER_BORDER_CLR);
                    end
                end
            end
            
            function clusterCallback(h)
                if h.getSelectedIndex<1
                    if ~isempty(dbm)
                        dbm.removeBorders;
                    end
                    dbm=[];
                else
                    dtls=Density.DETAILS(1:end-1);
                    detail=dtls{h.getSelectedIndex};
                    noBorders=isempty(hClusterBorders) ...
                        || ~ishandle(hClusterBorders(1));                    
                    if noBorders && ~isempty(dbm)
                        dbm.detail='';
                    end
                    if isempty(dbm) || ~isequal(dbm.detail, detail)
                        clusterIdentifiers=cluster(detail);
                        this.clusterIdentifiers=clusterIdentifiers;
                        showClusterLabels(4);
                    end
                end
                setClusterCount;
                showHeatMap(true);
            end
            
            function clues=cluster(detail, pu)
                if ~isempty(sprv) && isSupervising
                    newDetail=~isequal(sprv.clusterDetail, detail);
                    newSize=size(reduction, 1) ~= length(sprv.lastClusterIds);
                    noBorders=isempty(hClusterBorders) ...
                        || ~ishandle(hClusterBorders(1));                    
                    if ~newDetail && ~newSize && ~noBorders
                        clues=sprv.lastClusterIds;
                        return;
                    end
                    if ~newSize && ~isempty(legendRois)
                        legendRois.showBusyReclustering;
                    end
                    sprv.storeClusterIds=true;
                    if nargin<2
                        pu=PopUp('Changing the cluster detail level', 'south');
                    elseif islogical(pu)
                        pu=[];
                    else
                        pu.setText('Changing the cluster detail level');
                    end
                    if ~isempty(dbm)
                        dbm.removeBorders;
                    end
                    if newDetail || noBorders
                        sprv.setClusterDetail(detail)
                        if sprv.matchType>2 || ~beGraphic
                            sprv.findClusters(reduction, pu);
                        end
                        refreshUst(reduction, sprv.matchType, pu, false)
                        dbm=sprv.density;
                    elseif newSize
                        [~, ~, dbm]=sprv.findClusters(reduction, pu);
                    end
                    clues=sprv.lastClusterIds;
                    if beGraphic && args.n_components==2
                        hClusterBorders=dbm.drawBorders(curAxes, ...
                                    SuhRunUmap.CLUSTER_BORDER_CLR);
                    end
                    if nargin<2
                        pu.close;
                    end
                    if ~newSize && ~isempty(legendRois)
                        legendRois.hideBusy;
                    end                    
                elseif isempty(dbm) || ~isequal(dbm.detail, detail)
                    if nargin<2
                        pu=PopUp('Changing the cluster detail level', 'south');
                    elseif isa(pu, 'PopUp')
                        pu.setText('Changing the cluster detail level');
                    end
                    if ~isempty(dbm)
                        dbm.removeBorders;
                    end
                    [mins, maxs]=Supervisors.GetMinsMaxs(reduction);
                    if ~strcmpi('dbm', args.cluster_method_2D)
                        %using DBSCAN
                        [~, clues, dbm]=Density.FindClusters(reduction, ...
                            detail, args.cluster_method_2D, pu, ...
                            args.epsilon, args.minpts, ...
                            args.dbscan_distance, mins, maxs);
                    else
                        [~, clues, dbm]=Density.GetClusters(reduction, ...
                            detail, mins, maxs);
                    end
                    if beGraphic && args.n_components==2
                        hClusterBorders=dbm.drawBorders(curAxes, ...
                                    SuhRunUmap.CLUSTER_BORDER_CLR);
                    end
                    if nargin<2
                        pu.close;
                    end
                elseif nargout>0
                    if ~exist('pu', 'var')
                        pu=[];
                    end
                    [mins, maxs]=Supervisors.GetMinsMaxs(reduction);
                    [~, clues, dbm]=Density.FindClusters(reduction, ...
                        detail, args.cluster_method_2D, pu, ...
                        args.epsilon, args.minpts, ...
                        args.dbscan_distance, mins, maxs);
                    if beGraphic && args.n_components==2
                        hClusterBorders=dbm.drawBorders(curAxes, ...
                                    SuhRunUmap.CLUSTER_BORDER_CLR);
                    end
                    
                end
            end

            function showUstOptions(h)
                jMenu=PopUp.Menu;
                legendRois.addUstMenus(jMenu);
                jMenu.show(h, 5, 5);
            end
            
            function doLegendRoiButton
                prior=legendRois;
                if ~isempty(sprv)
                    if isempty(prior)
                        if ~isempty(btnUst)
                            btnUst.setVisible(isSupervising);
                        end
                    end
                    if ~isempty(sprv.plots)
                        legendRois=LabelLegendRoi(reduction, @roiMoved,...
                            sprv.plots.javaLegend, sprv.btns, ...
                            sprv.btnLbls, [], roiMap, args, this);
                        if isSupervising
                            if isempty(testSetLabels)
                                legendRois.setUst(...
                                    sprv, @refreshUst, inData);
                            else
                                legendRois.setUst(...
                                    sprv, @refreshUst, inData,...
                                    @refreshMatch);
                            end
                        end
                        if isempty(extras.supervisorMatchedLabels)
                            legendRois.labels=labels;
                        else
                            legendRois.labels=extras.supervisorMatchedLabels;
                        end
                        if ~argued.contains('roi_percent_closest')
                            if args.fast_approximation
                                legendRois.percentClosest=.91;
                            end
                        end
                    end
                end
                if ~isempty(legendRois)
                    legendRois.updateLegendGui(curAxes,...
                        args.roi_table~=0 && args.roi_table~=2, ... %create btn
                        fig, 'north east++', true);
                    legendRois.setSaveInfo(['umap_' args.reductionType...
                        '_roi' args.output_suffix '.properties' ], args.save_output,...
                        args.output_folder);
                    if ~isempty(prior)
                        legendRois.resetFromPrior(prior);
                    else
                        Gui.SetJavaVisible(legendRois.javaLegend);
                    end
                end
            end
            
            function refreshMatch(matchStrategy)
                locate_fig={fig, 'west++', true};
                if matchStrategy==2
                    word='F-measure';
                else
                    word='similarity';
                end
                set(0, 'CurrentFigure', fig);
                pu=PopUp(['Refreshing match based on ' word] );
                qft=sprv.qfDissimilarityTestSetPrior(...
                    reduction, inData, testSetLabels, ...
                    false, ... %sNOT with training
                    locate_fig, pu,...
                    [], matchStrategy, ...
                    false);
                if ~isempty(qft)
                    qft.setPredictionListener(@notifyPredictions);
                end
                pu.close;
            end

            function tryToRefreshSprvBtns(oldJavaBtns)
                try
                    if isempty(oldJavaBtns)
                        return;
                    end
                    it=oldJavaBtns.iterator;
                    newBtns=sprv.btns;
                    nms=sprv.plots.names;
                    nNames=length(nms);
                    seek=java.lang.String('&bull;</font></font>');
                    while it.hasNext
                        btn=it.next;
                        if ~btn.isSelected
                            bStr=btn.getText;
                            if isempty(bStr)
                                continue;
                            end
                            bStr=btn.getText;
                            %find portion of button name that does NOT
                            %change for reasons of size of subset
                            idx=bStr.indexOf('<b> training ');
                            if idx<5
                                %not likely to match ... OH well
                                %we DID try .... user will have to deselect
                                warning('No "<b> training"; in %s', ...
                                    char(bStr));
                                continue;
                            end
                            bStr=bStr.substring(0, idx);
                            nit=newBtns.iterator;
                            while nit.hasNext
                                newBtn=nit.next;
                                if newBtn.getText.startsWith(bStr)
                                    newBtn.setSelected(false);
                                    idx=bStr.indexOf(seek);
                                    if idx<0
                                        bStr=char(bStr);
                                        warning('No &bull; in %s', bStr);
                                        bStr=char(edu.stanford.facs.swing.Basics.RemoveXml(bStr));
                                    else
                                        bStr=char(bStr.substring(...
                                            idx+seek.length));
                                    end
                                    for bi=1:nNames
                                        if startsWith(nms{bi}, bStr)
                                            set(sprv.plots.Hs(bi), 'visible', 'off');
                                            break;
                                        end
                                    end
                                    break;
                                end
                            end
                        end
                    end
                catch ex
                    ex.getReport
                end
            end

            function [testSet, qfForClusterMatch]=refreshUst(...
                    data, matchType, pu, firstTime)
                if nargin<3
                    pu=[];
                end
                if ~firstTime
                    old=sprv.btns;
                end
                if beGraphic
                    [testSet, extras.supervisorMatchedLabels, qfForClusterMatch]=...
                        sprv.plotTestSet(umap, curAxes, data, pu, ...
                        matchType, true, true);
                    if ~firstTime
                        doLegendRoiButton;
                        tryToRefreshSprvBtns(old);
                    end
                end
            end

            function wbu
                wbUp=true;
                if ~seekingDataIsland
                    if ~RoiUtil.CanDoNew
                        if ~isempty(roi)
                            roiMoved;
                        end
                    end
                else
                    if ~isempty(hClusterBorders) ...
                        && ~ishandle(hClusterBorders(1))
                        hClusterBorders=dbm.drawBorders(curAxes, ...
                                    SuhRunUmap.CLUSTER_BORDER_CLR);
                    end
                    [roi_,clue]=dbm.getROI(curAxes, hClusterBorders, @roiMoved);
                    if ~isempty(dbmROI)
                        delete(dbmROI);
                    end
                    dbmROI=roi_;
                    if isempty(roi_)
                        if clue==0
                            msg('<html>No significant <br>cluster here ...</html>', ...
                                5, 'north east', 'Ooops...?', 'polygonGate.png');
                        end
                    else
                        RoiUtil.SetColor(roi, RoiUtil.EDIT_COLOR);
                        roi=dbmROI;
                        if ~isempty(jdDbm)
                            try
                                Gui.Shake(jdDbm);
                            catch
                            end
                        end
                    end
                end
            end
            
            function wbd
                if tb.jToolbar.isEnabled
                   if ~seekingDataIsland && ~drawingRoi
                       if ~isempty(hClusterBorders) ...
                               && ishandle(hClusterBorders(1)) ...
                               && strcmpi(get(hClusterBorders(1), ...
                                  'visible'), 'on')
                              cp=get(curAxes,'CurrentPoint');
                              cp=cp(1, 1:end-1);
                              nRoi=roiMap.size;
                              fprintf('# of ROIs is %d\n', nRoi);
                              if nRoi<2 %~any(RoiUtil.GetRows( roi, cp))
                                  if all(cp>=mins) && all(cp<=maxs)
                                      createRoiWithClusterPick(true);
                                      drawnow;
                                  end
                              end
                       end
                       if ~seekingDataIsland
                           %Gui.Shake(btnClusterRoi, 2, 'Pick clustering detail');
                           Gui.Shake(lblClu, 2, 'Click here to pick cluster polygons');
                       end
                   end
                end
                wbUp=false;
            end
            
            function qft=findQfMatch(detail)
                qft=[];
                qfds=this.extras.qfd;
                for i=1:length(qfds)
                    if strcmp(qfds{i}.contextDescription, detail)
                        qft=qfds{i};
                        return;
                    end
                end
            end

            function createRoiWithClusterPick(selectingCluster)
                if seekingDataIsland
                    jdDbm.setVisible(true);
                    setAlwaysOnTopTimer(jdDbm, 2, true);
                    return;
                end
                if isempty(dbm)
                    [~,~,dbm]=Density.FindClusters(reduction, 'medium');
                    setClusterDetail;
                end
                if (nargin==0 || ~selectingCluster) && ~isempty(args.save_roi)
                    unmatchedOnly=false;
                    if ~isempty(extras.qfd) && ~isempty(testSetLabels) ...
                            && strcmp(args.reductionType, UMAP.REDUCTION_BASIC)
                        qft=findQfMatch(dbm.detail);
                        if isempty(qft)
                            if askYesOrNo(sprintf(['<html>Run QFMatch for <b>%s</b>' ...
                                    ' level of cluster detail?</html>'], dbm.detail))
                                busy=Gui.ShowBusy(fig, Gui.YellowSmall('Re-running QFMatch'));
                                testBasicReduction(dbm.detail);
                                qft=findQfMatch(dbm.detail);
                                Gui.HideBusy(fig, busy);
                            end
                        end
                        if isempty(qft)
                            warning('QFMatch detail %s, current detail %s', args.cluster_detail{1}, dbm.detail)
                        else
                            N=length(qft.qf.sIds);
                            if N ~= dbm.numClusters
                                warning('%d QFMatch clusters and %d current clusters', N, dbm.numClusters);
                            else
                                sIds=qft.getUnmatched(true);
                                unmatchedOnly=~isempty(sIds);
                                if ~unmatchedOnly
                                    warning('No unmatched clusters!');
                                end
                            end
                        end
                    end
                    if unmatchedOnly
                        prop='SuhRunUmap.ConvertUnmatchedClusters';
                        choice=Gui.Ask(struct('where', 'north+', ...
                            'msg', 'What do you wish to do?'), ...
                            {'Convert all clusters to polygon gates.', ...
                            sprintf('Convert %d unmatched clusters', length(sIds)),...
                            'Let me select & convert individually'}, ...
                            prop, ...
                            'Clusters to polygons', 3);
                        if isempty(choice)
                            return;
                        end
                        if choice==1
                            convertAllClusters;
                            return
                        elseif choice==2
                            convertAllClusters(sIds);
                            return;
                        end

                    else
                        prop='SuhRunUmap.ConvertAllClusters';
                        choice=Gui.Ask(struct('where', 'north+', ...
                            'msg', 'What do you wish to do?'), ...
                            {'Convert all clusters to polygon gates.', ...
                            'Let me select & convert individually'}, ...
                            prop, ...
                            'Clusters to polygons', 2);
                        if isempty(choice)
                            return;
                        end
                        if choice==1
                            convertAllClusters;
                            return
                        end
                    end
                end
                seekingDataIsland=true;
                tb.setEnabled(false); 
                dbmROI=[];
                jdDbm=Density.MakePolgygonMsg(@(m,j)concludeDbm(m), fig);
                globals.showToolTip(btnClusterRoi, Html.WrapHr([...
                    'Click on any "data island" (cluster)'...
                    '<br>to create a polygon ROI!']), 10, 20);
                showClusterLabels(3);
            end

            function convertAllClusters(clusterIds)
                if isempty(args.save_roi)
                    return
                end
                clues=this.clusterIdentifiers;
                if nargin<1                    
                    clusterIds=unique(clues);
                    clusterIds(clusterIds==0)=[];
                end
                nClues=length(clusterIds);
                if beGraphic
                    txt=String.Pluralize2('cluster', nClues');
                    busy=Gui.ShowBusy(fig, Gui.YellowSmall(sprintf('Creating %d polygon gates', nClues)));
                    hWait=Gui.ProgressBar(sprintf('Converting %s to polygons', txt), true, fig, 'south++');
                end
                for i=1:nClues
                    if beGraphic
                        stop=getappdata(hWait,'canceling');
                        if stop
                            break;
                        end
                    end
                    clue=clusterIds(i);
                    clump=this.reduced_data(clues==clue,:);
                    [theRoi, ~]=RoiUtil.NewForXy(...
                        curAxes, clump, .00121, @roiMoved);
                    k=[clue '.' num2str(i)];
                    [n, named]=nameTheRegion(theRoi, false);
                    strClue=[ dbm.detail ' cluster #' num2str(i)  ];
                    if named
                        n=[n '( ' strClue  ')'];
                    else
                        n=strrep(n, args.reduction_algorithm, strClue);
                    end
                    if roiMap.containsKey(k)
                            fprintf('Cluster %s already gated...\n', n);
                        continue;
                    end
                    otherKey=LabelLegendRoi.FindName(roiMap, n);
                    if ~isempty(otherKey) && ~isequal(otherKey, k)
                        fprintf('Cluster %s already gated...\n', n);
                        continue;
                    end
                    roiMap.set(k, theRoi);
                    roiMap.set([k '.name'], n);
                    feval(args.save_roi, k, ...
                        theRoi, n, this.reduced_data, args);
                    if beGraphic
                        waitbar2a(i/nClues, hWait, ...
                            sprintf('%d of %s done', i, txt));
                    end
                end
                if beGraphic
                    Gui.HideBusy(fig, busy);
                    delete(hWait);
                end
            end
            
            function ok=concludeDbm(makePolygon)
                ok=false;
                try
                    Density.UnselectBorders(hClusterBorders);
                catch ex
                    ex.getReport
                end
                if ~makePolygon
                    if ~isempty(dbmROI)
                        delete(dbmROI);
                        roi=[];
                    end
                    seekingDataIsland=false;
                    ok=true;
                    tb.setEnabled(true); 
                elseif isempty(dbmROI)
                    nClues=length(unique(clusterIdentifiers));
                    if any(clusterIdentifiers==0)
                        nClues=nClues-1;
                    end
                    txt=sprintf([ ...
                            'No clusters are chosen ...<br><br>' ...
                            'Create %d polygons for ALL %d clusters?'],...
                            nClues, nClues);
                    if askYesOrNo(Html.WrapHr(txt), ...
                            'Convert ALL???', 'east+')
                        convertAllClusters;
                        return
                    end
                    return;
                end
                dbmROI=[];
                if ~isempty(args.save_roi)
                    if makePolygon
                        name=strtrim(nameTheRegion);
                    end
                    if makePolygon && isempty(name)
                        delete(roi);
                        roi=[];
                    elseif makePolygon
                        tb.setEnabled(false); %guard against re-entrancy ...
                        Gui.ShowFacs(fig, 'Saving region of interest');
                        saveNamedRegions(name);
                        RoiUtil.SetColor(roi, RoiUtil.NEW_COLOR);
                        Gui.HideBusy(fig);
                    end
                elseif ~isempty(roi)
                    RoiUtil.SetColor(roi, RoiUtil.NEW_COLOR);
                end
            end

            function createRoi(~, roiType)
                tb.setEnabled(false); 
                drawingRoi=true; %guard against re-entrancy ...
                roi=RoiUtil.New(curAxes, roiType, @roiMoved);
                drawingRoi=false; 
                tb.setEnabled(true);
                if ~isempty(args.save_roi)
                    name=strtrim(nameTheRegion);
                    if isempty(name)
                        msg(Html.WrapSmall([ ...
                            '<center>Deleted region of ' ...
                            'interest.<br>(Name was needed)' ...
                            '<hr></center>']), ...
                            4, 'north east+', 'Note...', ...
                            'umap.png');
                        delete(roi);
                        roi=[];
                    else
                        saveNamedRegions(name);
                    end
                else
                    txt=getTrainers(roi);
                    if ~isempty(txt)
                        html=Html.Wrap([txt '<hr>Note:  <i>This region '...
                            'does not have a name</i>'], 200);
                    else
                        html='';
                    end
                    if ~isempty(roi)
                        if ~isempty(html)
                            shakeNameBtn(true, html);
                        else
                            shakeNameBtn(true);
                        end
                    end
                end
            end

            function showRoiTable(rows, htmlRoiName)
                try
                    needToMake=isempty(roiTable) || ~ishandle(roiTable.table.table.fig);
                catch
                    needToMake=true;
                end
                if needToMake
                    if ~isempty(legendRois) %avoid covering up legend for labels
                        where='south east++';
                    else
                        where='east++';
                    end
                    roiTable=Kld.Table(inData(rows, :), parameter_names, ...
                        args.roi_scales,fig, htmlRoiName, 'south', ...
                        'Parameter', args.reduction_algorithm, ...
                        false, [], {fig, where, true}, false, densityBars);
                else
                    roiTable.refresh(inData(rows,:), htmlRoiName, rows);
                end
                figure(roiTable.getFigure);
                figure(fig);
            end
            
            function roiMoved(pRoi)
                if ~isempty(legendRois) && legendRois.buildingMultiplePolygons
                    return;
                end
                if nargin<1
                    pRoi=roi;
                else
                    roi=pRoi;
                end
                txt=getTrainers(pRoi);
                if ~isempty(txt)
                    shakeNameBtn(true, Html.Wrap(txt, 200));
                end
                try
                    if ~RoiUtil.IsNewRoi(pRoi)
                        if ~wbUp %necessary before r2018b
                            return;
                        end
                        newPos=RoiUtil.Position(pRoi);
                        if isequal(newPos, lastRoiPos)
                            return;
                        end
                        lastRoiPos=newPos;
                    end
                    try
                        rows=RoiUtil.GetRows(pRoi, reduction);
                    catch ex
                        ex.getReport
                        rows=[];
                    end
                    nRowsEnclosed=sum(rows);
                    fprintf('ROI found %d rows\n', nRowsEnclosed);
                    [~, roiName]=LabelLegendRoi.Find(roiMap, pRoi);
                    deleted=false;
                    if nRowsEnclosed<1
                        if isempty(roiName)
                            wording='region with no name';
                        else
                            wording=sprintf('region named "%s"', roiName);
                        end
                        if askYesOrNo(Html.WrapHr(sprintf('Delete this empty<br>%s?', ...
                                wording)), 'Region of interest??...', ...
                                'north+', true, '', 'run_umap.DeleteRoi')
                            LabelLegendRoi.Delete(roiMap, pRoi);
                            deleted=true;
                        end
                    end
                    if ~isempty(args.save_roi) && ~deleted
                        if ~isempty(roiName)
                            saveNamedRegions(roiName);
                        end
                    end
                    if ~cbSyncKld.isSelected
                        shakeSyncBtn;
                        return;
                    end

                    Gui.setEnabled(fig, false)
                    tb.setEnabled(false);
                    pu_=PopUp(Html.Wrap(['Synchronizing with '...
                        '<br>ParameterExplorer&nbsp;&nbsp;' ...
                        Html.ImgXy('pseudoBarHi.png', [], 1.2)]), ...
                        'north east+', 'One moment...', true, true);
                    showRoiTable(rows, String.ToHtml(roiName));
                    Gui.setEnabled(fig, true)
                    tb.setEnabled(true);
                    pu_.close;
                    drawnow
                    if pu_.cancelled
                        cbSyncKld.setSelected(false);
                        globals.showToolTip(cbSyncKld, Html.WrapHr(...
                            ['Synchronizing was <b>cancelled</b>...<br>'...
                            '(Click shaking button <b>above</b> to resume it)']), ...
                            -22, 23, 0, [], true, .31);
                        shakeSyncBtn;
                    else
                        if isempty(txt) && isempty(roiName) && isempty(args.save_roi)
                            MatBasics.RunLater(@(h,e)shakeNameBtn(true), 3);
                        end
                    end
                catch ex
                    ex.getReport
                    if exist('pu_', 'var')
                        pu_.close;
                        Gui.setEnabled(fig, true)
                        tb.setEnabled(true);
                    end
                end
            end
            
            function syncKld
                prop='run_umap.SyncKld';
                if cbSyncKld.isSelected
                    globals.set(prop, 'true');
                    if ~isempty(roi)
                        try
                            roiMoved
                        catch ex
                            ex.message
                        end
                    end
                else
                    globals.set(prop, 'false');
                end
                if ~isempty(roiTable)
                    roiFig=roiTable.getFigure;
                    if ishandle(roiFig)
                        figure(roiFig);
                        figure(fig);
                    end
                end
                
            end
            
            function testBasicReduction(details)
                if isempty(testSetLabels)
                    return
                end
                if ~strcmp(args.reductionType, UMAP.REDUCTION_BASIC)
                    return;
                end
                if nargin<1
                    details=args.cluster_detail;
                    createdPu=false;
                    if ~exist('pu', 'var')
                        if beQuiet
                            pu=[];
                        else
                            createdPu=true;
                            pu=PopUp('Matching results', 'west+', [], false, ...
                                [],[],false, args.parent_popUp);
                        end
                    end
                else
                    if ischar(details)
                        details={details};
                    end
                    createdPu=true;
                    pu=PopUp('Matching results', 'west+', [], false, ...
                        [],[],false, args.parent_popUp);
                end
                scenarios=args.match_scenarios;
                if nargin>0
                else
                    if args.match_predictions || args.confusion_chart
                        scenarios(end+1)=5;
                    end
                end
                nCluDtls=length(details);
                nScenarios=length(scenarios);
                last3or4=[];
                for c=1:nCluDtls
                    for ms=1:nScenarios
                        scenario=scenarios(ms);
                        reportProgress(sprintf(...
                            'Match clusters=%s, scenario=%d:"%s"', ...
                            details{c},scenario, ...
                            UmapUtil.GetMatchScenarioText(scenario, ...
                            args.reductionType)));
                        if scenario==3
                            matchStrategy=1; % match by QF dissimilarity
                        elseif scenario==4
                            matchStrategy=2; % match by F-measure overlap
                        elseif scenario==5
                            matchStrategy=2; % match by F-measure overlap
                        else
                            continue;
                        end
                        [clusterIds, numClusters, density]=UmapUtil.Cluster(...
                            reduction, details{c}, pu, ...
                            args.cluster_method_2D, args.minpts, ...
                            args.epsilon, args.dbscan_distance);
                        if strcmpi('dbm', args.cluster_method_2D)
                            dbm=density;
                        end
                        args.parameter_names=parameter_names;
                        if scenario==5
                            if ~isempty(pu)
                                pu.dlg.setTitle('Matching predictions to prior classification...');
                            end
                            if isempty(last3or4)
                                if beGraphic
                                    locate_fig={fig, 'north west++', true};
                                else
                                    locate_fig=false;
                                end
                                qft=UmapUtil.Match(args, inData, ...
                                    testSetLabels, labelMap, clusterIds,  ...
                                    details{c}, 1,...
                                    locate_fig, pu, true);
                            else
                                if args.match_predictions
                                    if beGraphic
                                        locate_fig={last3or4.fig, 'east+', true};
                                    else
                                        locate_fig=false;
                                    end
                                    pred=last3or4.getPredictionsOfThese;
                                    qft=pred.showTable(locate_fig, pu);
                                    qft.context=last3or4.context;
                                    qft.context.matchScenario=5;
                                end
                                if args.mds
                                    last3or4.mds;
                                end
                                if args.confusion_chart
                                    last3or4.confusionChart;
                                end
                                if args.false_positive_negative_plot
                                    last3or4.browseFalsePosNeg;
                                end
                                if ~isempty(args.match_webpage_file)
                                    last3or4.addToWebPage(...
                                        args.match_webpage_file)
                                end
                            end
                            if ~isempty(qft)
                                qft.qf.description='UMAP unsupervised predictions';
                                extras.qfd{end+1}=qft;
                            end
                        else
                            if ~isempty(pu)
                                pu.dlg.setTitle('Matching clusters to prior classification...');
                            end
                            [qft, tNames]=UmapUtil.Match(args, inData, ...
                                testSetLabels, labelMap, clusterIds,  ...
                                details{c}, matchStrategy, ...
                                false, pu, false, [], probability_bins);
                            if ~isempty(qft)
                                if scenario==3 || scenario==4
                                    last3or4=qft;
                                end
                                qft.qf.description='UMAP unsupervised';
                                extras.qfd{end+1}=qft;
                                if beGraphic
                                    if args.match_table_fig
                                        if ~isempty(qft.fig)
                                            SuhWindow.Follow(qft.fig, fig, ...
                                                'north west++', true);
                                            SuhWindow.SetFigVisible(qft.fig);
                                        end
                                    end
                                    if args.match_histogram_figs
                                        if ~isempty(qft.qHistFig)
                                            SuhWindow.SetFigVisible(qft.qHistFig);
                                        end
                                        if ~isempty(qft.fHistFig)
                                            SuhWindow.SetFigVisible(qft.fHistFig);
                                        end
                                    end
                                    matchedLbls=UmapUtil.GetMatches(reduction, qft.qf, ...
                                        tNames, labelMap, density, clusterIds, numClusters);
                                    [plots, qft.btns, qft.btnLbls]...
                                        =Supervisors.Plot( reduction, matchedLbls,...
                                        labelMap, nCols-labelCols, umap, curAxes, ...
                                        true, false, true, args);
                                    legendRois=LabelLegendRoi(...
                                        reduction, @roiMoved, ...
                                        plots.javaLegend, qft.btns, ...
                                        qft.btnLbls, matchedLbls, ...
                                        roiMap, args, this);
                                    if ~argued.contains('roi_percent_closest')
                                        if args.fast_approximation
                                            legendRois.percentClosest=.91;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                if ~isempty(extras.qfd)
                    extras.doMatchOutput(nCols-labelCols);
                    if beGraphic
                        extras.seeMatches(args.match_html);
                    end
                end
                if createdPu
                    pu.close;
                end
            end
            
            
            function [map, halt]=getLabelMap(lbls)
                halt=false;
                map=[];
                if isempty(args.label_file)
                    warning(['label_column without label_file '...
                        'to match/supervise, will use default names/colors']);
                    args.buildLabelMap=true;
                end
                if args.buildLabelMap
                else
                    if exist(args.label_file, 'file')
                        map=File.ReadProperties( args.label_file);
                        if isempty(map)
                            problem='load';
                        end
                    elseif ~isempty(args.label_file)
                        problem='find';
                    end
                    if isempty(map)
                        if askYesOrNo(['<html>Cannot ' problem ' the '...
                                ' label file <br><br>"<b>' globals.smallStart ...
                                args.label_file globals.smallEnd '</b>"<br><br>'...
                                '<center>Use default names & colors?</center>'...
                                '<hr></html>'], 'Error', 'north west', true)
                            args.buildLabelMap=true;
                        else
                            halt=true;
                        end
                    end
                end
                if args.buildLabelMap
                    map=java.util.Properties;
                    u=unique(lbls)';
                    nU=length(u);
                    if nU/nRows > .2
                        if ~acceptTooManyLabels(nU)
                            halt=true;
                            return;
                        end
                    end
                    for i=1:nU
                        key=num2str(u(i));
                        map.put(java.lang.String(key), ['Subset #' key]);
                        map.put([key '.color'], num2str(Gui.HslColor(i, nU)));
                    end
                end
                if args.color_defaults
                    ColorsByName.Override(map, args.color_file, beQuiet);
                end
            end
            
            function clues=doClusters()
                if isempty(reduction)
                    clues=[];
                else
                    if ~beGraphic
                        pu=false;
                    end
                    clues=cluster(args.cluster_detail{1}, pu);
                    if strcmpi(args.cluster_output, 'graphic')
                        if ~isempty(clues)
                            if ~exist('xLabel', 'var')
                                dimInfo=sprintf('  %dD\\rightarrow%dD', nCols-labelCols, ...
                                    args.n_components);
                                xLabel=[args.reduction_algorithm '-X' dimInfo];
                                yLabel=[args.reduction_algorithm '-Y' dimInfo];
                                zLabel=[args.reduction_algorithm '-Z' dimInfo];
                            end
                            cp=ClusterPlots.Go(fig, reduction, clues, [], xLabel, ...
                                yLabel, zLabel, true, [], false, true, false,...
                                true, 'south west++', args);
                            if args.n_components==2
                                if isequal('dbscan', args.cluster_method_2D)
                                    clue='dbscan';
                                else
                                    clue='dbm';
                                end
                            else
                                clue='dbscan';
                            end
                            annotateClues(get( ...
                                cp.ax, 'Parent'), ...
                                args.cluster_detail{1}, ...
                                clue, args.epsilon, args.minpts, ...
                                args.dbscan_distance);
                        end
                    end
                end
            end
            
            function lbl=annotateClues(fig, detail, clue, epsilon, minpts, dist)
                X=.005;
                Y=.872;
                W=.58;
                H=.115;
                [epsilon, minpts]=Density.GetDbscanParameters(detail, ...
                    epsilon, minpts);
                info=['clue method="', clue '", detail="' detail '"'];
                info2=['epsilon=' num2str(epsilon) ', minpts=' num2str(minpts) ...
                    ', DBSCAN distance="' dist '"'];
                lbl=Gui.TextBox({['\color{blue} ' info], ...
                    ['\fontsize{10} ' info2]}, fig, ...
                    'position', [X Y W H],...
                    'fontSize', 11);
            end
            
            function hContour=updatePlot(data, lastCall, ...
                    contours, pseudoColor, contourPercent, reset)
                if nargin<6
                    reset=true;
                    if nargin<5
                        contourPercent=10;
                        if nargin<4
                            contours=true;
                            if nargin<3
                                pseudoColor=true;
                            end
                        end
                    else
                        args.contour_percent=contourPercent;
                    end
                end
                hContour=[];
                doingDensity3D=false;
                labelsDone=true;
                if nargin<2
                    lastCall=false;
                end
                if ~lastCall && nargin<4
                    contours=false;
                end
                if nargin>0
                    if isempty(xLabel)
                        dimInfo=sprintf('  %dD\\rightarrow%dD', nCols-labelCols, ...
                            args.n_components);
                        xLabel=[args.reduction_algorithm '-X' dimInfo];
                        yLabel=[args.reduction_algorithm '-Y' dimInfo];
                        if args.n_components>2
                            zLabel=[args.reduction_algorithm '-Z' dimInfo];
                        end
                    end
                    if args.n_components>2
                        nD=size(data, 2);
                        assert(nD==args.n_components);
                        if ~plotLabels(data, lastCall)
                            labelsDone=false;
                            if args.frequencyDensity3D && pseudoColor
                                doingDensity3D=true;
                                Gui.PlotDensity3D(curAxes, data, 64, 'iso',...
                                    xLabel, yLabel, zLabel, args);
                            else
                                Gui.PlotNeighDist3D(curAxes, data, ...
                                    args.n_neighbors);
                            end
                        end
                        if args.n_components>3
                            title(curAxes, ['NOTE:  Only 3 of \color{red}' ...
                                num2str(args.n_components) ...
                                ' dimensions being shown...']);
                        end
                    else
                        if pseudoColor && plotLabels(data, lastCall)
                            if ~isempty(paramAnnotation)
                                if ishandle(paramAnnotation)
                                    set(paramAnnotation, 'visible', 'off');
                                end
                            end
                        else
                            labelsDone=false;
                            if lastCall
                                hContour=ProbabilityDensity2.Draw(curAxes, data, ...
                                    contours, pseudoColor, reset, .05, ...
                                    contourPercent, args);
                            else
                                hContour=ProbabilityDensity2.Draw(curAxes, data, ...
                                    contours, pseudoColor, reset, 0, ...
                                    contourPercent,args);
                            end
                        end
                    end
                end
                if ~labelsDone
                    if nargin>0
                        umap.adjustLims(curAxes, data );
                    else
                        umap.adjustLims(curAxes);
                    end
                    xlabel(curAxes, xLabel);
                    ylabel(curAxes, yLabel);
                    if args.n_components>2
                        zlabel(curAxes, zLabel);
                    end
                    grid(curAxes, 'on')
                    set(curAxes, 'plotboxaspectratio', [1 1 1])
                    if lastCall
                        if ~doingDensity3D
                            Gui.StretchLims(curAxes, data, .04);
                        end
                    end
                end
                if lastCall
                    if args.n_components==2
                        if isSupervising
                            if ~hasLabels
                                hClusterBorders=...
                                    sprv.drawClusterBorders(curAxes);
                            end
                        end
                    end
                end
                drawnow;
            end

            function ok=rePlotLabels()
                if ~isempty(matchedLbls)
                    Supervisors.Plot(reduction, matchedLbls, labelMap, ...
                        nCols-labelCols, umap, curAxes, true, ...
                        false, true, args);
                elseif hasLabels
                    ok=true;
                    Supervisors.Plot(reduction, labels, labelMap, ...
                        nCols-labelCols, umap, curAxes, true, ...
                        false, true, args);
                end
            end

            function ok=plotLabels(data, lastCall)
                if hasLabels
                    ok=true;
                    if lastCall
                        [plots, btns, btnLbls]=Supervisors.Plot( ...
                            data, labels, labelMap, nCols-labelCols, ...
                            umap, curAxes, true, false, true, args);
                        if ~isempty(plots) && ~isempty(umap.supervisors)
                            umap.supervisors.plots=plots;
                            umap.supervisors.btns=btns;
                            umap.supervisors.btnLbls=btnLbls;
                        end
                        Gui.StretchLims(curAxes, data, .04);
                        if ~isempty(sprv)
                            sprv.prepareForTemplate;
                            if args.qf_tree
                                sprv.inputData=umap.raw_data;
                                [~,qft]=sprv.qfTreeSupervisors(...
                                    {fig, 'north west++', true}, ...
                                    [], 'UMAP training set');
                                if ~isempty(qft) && ~isempty(qft.fig)
                                    extras.qft=qft;
                                end
                            end
                        end
                    else
                        Supervisors.Plot(data, labels, labelMap, nCols-labelCols, ...
                            umap, curAxes, false, false, true, args);
                    end
                elseif isSupervising
                    ok=true;
                    if lastCall
                        if isempty(sprv.embedding)
                            sprv.embedding=umap.embedding;
                        end
                        if firstPlot && args.see_training
                            [curAxes, ~,~,extras.supervisorMatchedLabels]...
                                =sprv.plotTrainingAndTestSets(...
                                data, curAxes, umap, pu, args.match_supervisors(1), ...
                                true);
                        else
                            [~,firstQf]=...
                                refreshUst(data, args.match_supervisors(1), pu, true);
                        end
                        doQfs(data);
                        drawnow;
                    else
                        if firstPlot && args.see_training
                            curAxes=sprv.plotTrainingAndTestSets(...
                                data, curAxes, umap, pu, progressMatchType, false);
                        else
                            sprv.plotTestSet(umap, curAxes, ...
                                data, pu, progressMatchType, false);
                        end
                        firstPlot=false;
                    end
                else
                    ok=false;
                end
            end
            
            function doQfs(data)
                scenarios=args.match_scenarios;
                if args.match_predictions || args.confusion_chart
                    if isequal(scenarios,0)
                        scenarios=5;
                    else
                        scenarios(end+1)=5;
                    end
                end
                nScenarios=length(scenarios);

                if args.qf_tree || all(scenarios>0)
                    creatingPu=~exist('pu', 'var') || isempty(pu);
                    if creatingPu
                        if beQuiet
                            pu=[];
                        else
                            pu=PopUp('Matching results', 'north west+', [], false);
                        end
                    end
                    sprv.inputData=umap.raw_data;
                    cascading={};
                    hasFig=exist('fig', 'var') && Gui.IsVisible(fig);
                    if hasFig
                        scrFig=fig;
                    else
                        scrFig=[];
                    end
                    if args.qf_tree
                        if ~beQuiet
                            disp('Computing QF-tree(s)');
                        end
                        if beGraphic
                            pu.dlg.setTitle('Creating QF-Tree...');
                            [~,qft]=sprv.qfTreeSupervisors(...
                                {fig, 'north west++', true}, pu, ...
                                'UMAP template''s training set');
                        else
                            [~,qft]=sprv.qfTreeSupervisors(false, pu, ...
                                'UMAP template''s training set');
                        end
                        if ~isempty(qft) && ~isempty(qft.fig)
                            extras.qftSupervisors=qft;
                            if beGraphic
                                [~,qft]=sprv.qfTreeSupervisees(data, ...
                                    inData, {fig, 'west++', true}, pu);
                            else
                                [~,qft]=sprv.qfTreeSupervisees(data, ...
                                    inData, false, pu);
                            end
                            if ~isempty(qft)
                                extras.qft=qft;
                            end
                        end
                    end
                    if all(scenarios>0)
                        if ~beQuiet
                            matchProgress('Matching UST results');
                        end
                        if hasFig
                            figure(scrFig);
                        end
                        matchTypes=unique(args.match_supervisors);
                        matchTypes=[args.match_supervisors(1) ...
                            matchTypes(matchTypes ~= args.match_supervisors(1))];
                        nMatchTypes=length(matchTypes);
                        nCluDtls=length(args.cluster_detail);
                        for c=1:nCluDtls
                            for mi=1:nMatchTypes
                                if c>1 || (mi>1 || ~beGraphic)
                                    matchType=matchTypes(mi);
                                    if c>1
                                        if matchType>=3 % nearest neighbor no clustering
                                            continue;
                                        elseif mi==1
                                            sprv.initClustering(...
                                                args.cluster_detail{c}, ...
                                                args.cluster_method_2D, ...
                                                args.minpts, ...
                                                args.epsilon, ...
                                                args.dbscan_distance);
                                            sprv.computeAndMatchClusters(data,...
                                                matchType, pu);
                                        end
                                    end
                                    if ~beQuiet
                                        if matchType<3
                                            word=[' @ "' args.cluster_detail{c} '"'];
                                        else
                                            word='';
                                        end
                                        matchProgress(...
                                            sprintf('New match type %s %s', ...
                                            UmapUtil.GetMatchTypeLongText(matchType, ...
                                            args.reductionType, args.n_components, ...
                                            nCols-labelCols), word), ...
                                            sprintf('Match type=%d, clu=%s', ...
                                            matchType, args.cluster_detail{c}));
                                    end
                                    if matchType==4
                                        sprv.changeMatchType(inData, ...
                                            matchType, pu);
                                    else
                                        sprv.changeMatchType(data, ...
                                            matchType, pu);
                                    end
                                end
                                last3or4=[];
                                for msi=1:nScenarios
                                    predictions=false;
                                    scenario=scenarios(msi);
                                    matchStrategy=1;
                                    if beGraphic
                                        locate_fig={fig, 'west++', true};
                                    else
                                        locate_fig=false;
                                    end                                    
                                    report=sprintf(...
                                        'Match scenario=%d:"%s"', ...
                                        scenario, UmapUtil.GetMatchScenarioText(...
                                        scenario, args.reductionType));
                                    if scenario==2
                                        if ~isempty(pu)
                                            pu.dlg.setTitle('Matching UMAP to training set...');
                                        end
                                        %match training set prior classification
                                        %to ust trained re-classification of test set
                                        [~,qfd]=sprv.qfDissimilarity(data, ...
                                            inData, locate_fig, pu, firstQf);
                                        
                                        firstQf=[];
                                        reportProgress(report);
                                    else
                                        if scenario==1
                                            %match training set classification to
                                            %prior classification of test set
                                            %ONLY needed once
                                            if mi>1 || c>1
                                                continue;
                                            end
                                            withTraining=true;
                                            if ~isempty(pu)
                                                pu.dlg.setTitle('Matching training/test set...');
                                            end
                                        else
                                            %match ust re-classification of test
                                            %set to prior classification of test set
                                            if ~isempty(pu)
                                                pu.dlg.setTitle('Matching UMAP to prior classification...');
                                            end

                                            withTraining=false;
                                            if scenario==4
                                                matchStrategy=2;
                                            elseif scenario==5
                                                matchStrategy=1;
                                                predictions=true;
                                                if ~isempty(pu)
                                                    pu.dlg.setTitle('Matching predictions to prior classification...');
                                                end
                                                if ~isempty(last3or4)
                                                    if args.match_predictions
                                                        if beGraphic
                                                            locate_fig={...
                                                                last3or4.fig,...
                                                                'north east+', true};
                                                        else
                                                            locate_fig=false;
                                                        end
                                                        pred=last3or4.getPredictionsOfThese;
                                                        qfd=pred.showTable(locate_fig, pu);
                                                        qfd.context=last3or4.context;
                                                        qfd.context.matchScenario=5;
                                                        if ~isempty(qfd)
                                                            extras.qfd{end+1}=qfd;
                                                        end
                                                    end
                                                    if args.mds
                                                        last3or4.mds;                                                        
                                                    end
                                                    if args.confusion_chart
                                                        last3or4.confusionChart;                                                        
                                                    end
                                                    if args.false_positive_negative_plot
                                                        last3or4.browseFalsePosNeg;
                                                    end
                                                    if ~isempty(args.match_webpage_file)
                                                        last3or4.addToWebPage(...
                                                            args.match_webpage_file)
                                                    end
                                                    continue;
                                                end
                                            else
                                                predictions=false;
                                            end
                                        end
                                        if isempty(testSetLabels)
                                            warning(...
                                                ['Cannot do qf dissimilarity if'...
                                                ' label_column is not provided']);
                                            continue;
                                        end
                                        reportProgress(report);
                                        qfd=sprv.qfDissimilarityTestSetPrior(...
                                            data, inData, testSetLabels, ...
                                            withTraining, locate_fig, pu,...
                                            [], matchStrategy, ...
                                            predictions);
                                        if scenario==3 || scenario==4
                                            last3or4=qfd;
                                        end
                                    end
                                    if ~isempty(qfd)
                                        extras.qfd{end+1}=qfd;
                                    else
                                        if creatingPu
                                            pu.close;
                                        end
                                        return;
                                    end                                    
                                end
                            end
                        end
                        extras.doMatchOutput(nCols-labelCols);
                        if beGraphic
                            if ~isempty(cascading)
                                Gui.CascadeFigs(cascading, false, true, 70, 2, ...
                                    true, false, scrFig, args.cascade_x);
                            end
                            if inputDataIsCsv && ischar(args.template_file)
                                h1=['<h3>' args.csv_file_or_data '<br>'...
                                    args.template_file '</h3>'];
                            elseif ischar(args.template_file)
                                h1=['<h3>' args.template_file '</h3>'];
                            else
                                h1=[];
                            end
                            if ~isempty(args.match_file)
                                extras.saveMatchFiles(h1);
                                if args.match_html==1
                                    extras.seeMatches(2, h1);
                                else
                                    extras.seeMatches(-1, h1);
                                end
                            else
                                extras.seeMatches(args.match_html, h1);
                            end
                        end
                        if hasFig
                            figure(scrFig);
                        end
                    end
                    if creatingPu
                        if ~isempty(pu)
                            pu.close;
                        end
                    end
                end
            end
            
            function matchProgress(s, ttl)
                if ~beQuiet
                    if isempty(args.parent_popUp) && exist('pu', 'var') ...
                            && ~isempty(pu)
                        if nargin>1
                            pu.dlg.setTitle(ttl);
                        else
                            pu.dlg.setTitle(s);
                        end
                    end
                    fprintf('%s %s\n',args.parent_context,  s);
                end
            end
            
            
            function keepComputing=progress_report(objectOrString)
                try
                    if ischar(objectOrString)
                        if ~ishandle(fig)
                            msg(Html.WrapHr(['Terminating since window has '...
                                'been closed ...<br>Avoid this by '...
                                'setting parameter <b><i>verbose</i></b> '...
                                'to ''<b>none</b>'' or ''<b>text</b>''!']), ...
                                8, 'south west');
                            keepComputing=false;
                            return;
                        end
                        if ~isequal(objectOrString, ...
                                StochasticGradientDescent.FINDING_ISLANDS)
                            drawnow;
                            if ~String.StartsWith(objectOrString, KnnFind.PROGRESS_PREFIX)
                                pu.pb.setValue(pu.pb.getValue+1);
                            end
                        end
                        keepComputing=~pu.cancelled;
                        pu.pb.setString(objectOrString);
                        pu.pack;
                        pu.showTimeSpent;
                        Gui.HideBusy(fig, busy);
                        simplicialSize=simplicialSize+.025;
                        switch simplicialIdx
                            case 1
                                suffix2='<br><br><br><br>';
                            case 2
                                suffix2='<br><br><br><br>';
                            otherwise
                                suffix2='<br><br><br><br>';
                        end
                        simplicialIdx=simplicialIdx+1;
                        busy=Gui.ShowBusy(fig, [Html.WrapSmallTags(...
                            char(objectOrString)) suffix2], ...
                            'simplicialComplex.png', simplicialSize);
                        return;
                    end
                    if objectOrString.getEpochsDone<10
                        Gui.HideBusy(fig, busy);
                    end
                    keepComputing=~pu.cancelled;
                    done=objectOrString.getEpochsDone-1;
                    toDo=objectOrString.getEpochsToDo;
                    pu.pb.setValue(3+(pu.pb.getMaximum*(done/toDo)));
                    pu.pb.setMaximum(toDo);
                    pu.pb.setString(sprintf('%d/%d epochs done', done, toDo));
                    if isvalid(runAnnotation)
                        delete(runAnnotation);
                    end
                    if args.rescale~=1
                        updatePlot(objectOrString.getEmbedding...
                            * args.rescale);
                        axis(curAxes, 'square');
                        drawnow;
                    else
                        updatePlot(objectOrString.getEmbedding);
                    end
                    pu.showTimeSpent;
                    return;
                catch ex
                    if umap.doing_stochastic_gradient_descent
                        throw(ex);
                    else
                        ex.getReport
                        keepComputing=false;
                    end
                end
            end
            
            
            function csvFile=downloadCsv
                csvFile=[];
                
                zip=fullfile(UmapUtil.LocalSamplesFolder, 'samples.zip');
                if ~isempty(WebDownload.GetZipIfMissing(zip, ...
                        WebDownload.ResolveUrl))
                    csvFile='sample10k.csv';
                    msg(Html.WrapHr(['Samples files are stored in<br>'...
                        UmapUtil.LocalSamplesFolder]), 8, 'south east+');
                end
            end
            
            
            function reportProgress(report, starting)
                if ~isempty(args.parent_popUp)
                    if nargin>1 && starting
                        args.parent_popUp.setText(['Test ' ...
                            args.parent_context ' ' args.description]);
                        args.parent_popUp.dlg.pack;
                    else
                        disp([num2str([args.parent_popUp.pb.getValue+1 args.parent_popUp.pb.getMaximum]) ' ' report]);
                        args.parent_popUp.incrementProgress;
                    end
                    args.parent_popUp.setText2(report);
                end
                if ~beQuiet
                    fprintf('%s %s\n',args.parent_context, report);
                end
                
            end
            
            function dispNoDbScan
                warning(['DBSCAN for clustering in 3+D is not available ... '...
                    '\nDownload from MathWorks File Exchange: '...
                    'https://www.mathworks.com/matlabcentral/fileexchange/52905-dbscan-clustering-algorithm']);
            end
            
            function ok=acceptTooManyLabels(nU)
                ok=true;
                txt=sprintf(['You have %s unique labels...<br>'...
                    'This is %s of the actual data rows...'], ...
                    String.encodeInteger(nU), String.encodePercent(...
                    nU/nRows, 1, 1));
                if ~beGraphic
                    ok=false; %#ok<NASGU>
                    error('This many labels is not supported, %s', txt);
                    
                end
                if ~askYesOrNo(Html.WrapHr(...
                        sprintf(['Interesting ....%s'...
                        '<br><br>So this then will be very SLOW ...'...
                        '<br><br><b>Continue</b>????'], ...
                        txt)))
                    ok=false;
                    return;
                end
            end
        end
        
    end

end