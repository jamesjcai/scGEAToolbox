classdef ArgumentClinic < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab 
%   License: BSD 3 clause

    properties(Constant)
        VERSION='4.09';
        PROP_UPDATE_CHECK_TIME='suh_pipelines.CheckUpdateTime';
    end
    
    properties(SetAccess=private)
        dlg;
        app;
        H;
        fjbFig;
        argsObjs=cell(1,3);
        jDlgClinics=cell(1,3);
    end
    
    methods
        function this=ArgumentClinic()
            persistent singleton;
            if ~isempty(singleton) %&& singleton.dlg.isVisible
                this=singleton;
                this.fjbFig=ArgumentClinic.FlowJoBridgeFig;
                this.dlg.setVisible(true);
                this.dlg.requestFocus;
                return;
            end           
            this.app=BasicMap.Global;
            BasicMap.SetHighDef(this.app, ...
                Gui.hasPcHighDefinitionAProblem);
            MatBasics.WarningsOff;
            pnl=ArgumentClinic.TopPanel(this.app,...
                @(idx)selectClinic(this, idx));
            this.fjbFig=ArgumentClinic.FlowJoBridgeFig;
                where='center';
            this.dlg=msg(pnl, 0, where,...
                [' Stanford University''s '...
                'Herzenberg pipelines (v' ...
                ArgumentClinic.VERSION ')']);
            dlg2=handle(this.dlg, 'CallbackProperties');
            set(dlg2, 'WindowClosingCallback', @(h,e)close());
            ArgumentClinic.CheckForUpdate(false);
            SuhWindow.Follow(this.fjbFig, this.dlg, ...
                'south west+', false);
            SuhWindow.SetFigVisible(this.fjbFig)
            Gui.ChimeUp;
            singleton=this;            
        end
        
        function selectClinic(this, idx)
            if ~isempty(this.jDlgClinics{idx})
                this.jDlgClinics{idx}.setVisible(true);
                this.jDlgClinics{idx}.requestFocus;
            elseif ~isempty(this.argsObjs{idx})
                argsObj=this.argsObjs{idx};
            else
                try
                    if idx==3
                        argsObj=SuhMatch.GetArgsWithMetaInfo();
                        where='south++';
                        forRunningWhat='for running QFMatch/QF-tree';
                        ttl='subset characterization';
                    elseif idx==1
                        [file, example]=SuhEpp.EliverArgs;
                        argsObj=SuhEpp.GetArgsWithMetaInfo(file, example{:});
                        where='south east++';
                        forRunningWhat='for running EPP';
                        ttl='unsupervised subset identification';
                    else
                        example={'label_column', 'end', ...
                            'match_scenarios', 2};
                        argsObj=UmapUtil.GetArgsWithMetaInfo(...
                            'eliverLabeled.csv',  example{:});
                        forRunningWhat='for running UMAP/UST';
                        ttl='unsupervised & supervised subset identification';
                        where='south west++';
                    end
                    
                    jd=javax.swing.JDialog;
                    jd.getContentPane.add(Gui.FlowPanelCenter(15,11,...
                        argsObj.getArgumentClinicPanel(forRunningWhat)));
                    jd.setTitle(['SUH ' ttl]);
                    jd.pack;
                    argsObj.javaWindow=jd;
                    this.jDlgClinics{idx}=jd;
                    SuhWindow.Follow(this.jDlgClinics{idx}, ...
                        this.dlg, where, true);
                    jd.setVisible(true)
                catch ex
                    msgError(['<html>' Html.Exception(ex, BasicMap.Global) '</html>']);
                    ex.getReport
                end
            end
        end
    end
    
    methods(Static)
        function tip=FlowJoStory
            tip=['<html><center>'...
                Html.WrapTable([ '<b>The FlowJoBridge!</b><'...
                '<hr>This integrates FlowJo workspaces with the ' ...
                'core algorithms of AutoGate.<br><br>' ...
                'Born at the Herzenberg Lab '...
                'FlowJo TM was primarily written by Mario Roederer '...
                'and Adam Treister, based on concepts developed '...
                'with David Parks, Martin Bigos, and Wayne Moore.'...
                '<br><br>FlowJo research & development is now led '...
                'by <b>Josef Spidlen</b> at <u>BD Life Sciences</u>.<br><br>'...
                'See <font color=''blue''>https://www.flowjo.com/about/company/story</font>'], 6, 200) ...
                '</center></html>'];
        end
        
        function pnl=TopPanel(app, fnc)
            if nargin<2
                app=BasicMap.Global;
                if nargin<1
                    fnc=[];
                end
            end
            [~,leonard]=Gui.ImageLabel(Html.Wrap(Html.ImgXy(...
                'Leonard.png', [], .5, false, false, app)), [], ...
                'See about our lab''s founder', @seeLen);
            [~,arthur]=Gui.ImageLabel(['<html>' ...
                '<center>Try our <u>new bridge</u> to <br>' ...
                Html.ImgXy('flowJoTM.png',[], .394, false,false, app) ...
                 '<br>' Html.WrapSmall('(<i>from BD Life Sciences</i>)')...
                '</html>'], [], ArgumentClinic.FlowJoStory, ...
                @(h,e)flowJoBridge(h));
            btnUpdate=Gui.NewBtn(Html.WrapSmallBold(...
            '<font color="blue">Check for<br>Updates</font>'),...
                @(h,e)ArgumentClinic.CheckForUpdate(true),...
                'Click to checker for newer version of pipelines');
            pipeline=Gui.Label(['<html>'...
                Html.ImgXy('pipeline.png',[], .94, false,false, app) ...
                '</html>']);
            wayne=html('wayneMoore2.png', .39, ...
                '<u>unsupervised</u> identification', ...
                'exhaustive projection pursuit (EPP)', 'Wayne Moore');
            connor=html('connor.png',.15, ...
                ['unsupervised & <u>supervised</u> '...
                '<br>&nbsp;&nbsp;&nbsp;identification'], ...
                'parameter reduction (UMAP/UST)', 'Connor Meehan');
            orlova= html('darya.png',.2, ...
                'characterization', 'QFMatch/QF-tree', 'Darya Orlova');
            items={ wayne, connor, orlova, connor};
            jRadio=Radio.PanelAndCallback(@resolve, true, items);
            normalFont=javax.swing.UIManager.getFont('Label.font');
            font=java.awt.Font(normalFont.getName, ...
                normalFont.ITALIC, normalFont.getSize+3);
            pnlLen=Gui.BorderPanel([], 0, 8, ...
                'North', leonard, 'South', ...
                Html.WrapSmallBold('Len Herzenberg'));
            pnlSouthWest=Gui.BorderPanel;
            pnlSouthWest.add(pnlLen, 'West');
            pnlSouthWestEast=Gui.BorderPanel([],0,4, 'North', ...
                Gui.Panel(arthur), 'South', Gui.Panel(btnUpdate));
            pnlSouthWest.add(pnlSouthWestEast, 'East');
            if app.highDef
                pnl=Gui.Panel( Gui.BorderPanel([], 2, 15, ...
                    'North', pipeline, 'South', ...
                    pnlSouthWest), Gui.BorderPanel([], 0, 2, 'North',...
                    Html.WrapHr(['<font color="blue"><b>'...
                    'Choose a data subsetting pipeline</b></font>']),...
                    'Center', jRadio));
            else
                pnl=Gui.Panel( Gui.BorderPanel([], 2, 15, ...
                    'North', pipeline, 'South', ...
                    pnlSouthWest), jRadio);
                Gui.SetTitledBorder('Choose a data subsetting pipeline', ...
                    pnl, font, java.awt.Color(.3, .3, .7), 2);
            end
            function str=html(img, scale, words, via, provider)
                img=Html.ImgXy(img,[], scale, ...
                    false,false,app);
                str=Html.Wrap(['<table cellspacing="5"><tr><td>' ...
                    img '</td><td><font color="blue"><b>Subset '...
                    words '</b></font><br>via <i>', via ...
                    '</i><br>' Html.WrapBoldSmall([' ' provider])...
                    '</td></tr></table>']);
            end
            
            function flowJoBridge(h)
                FlowJoTree.Open(Gui.WindowAncestor(h));
            end
            
            function seeLen(h,e)
                web('https://www.pnas.org/content/110/52/20848', '-browser');
            end
            
            function resolve(h,e)
                item=char(e.getActionCommand);
                idx=StringArray.IndexOf(items, item);
                if isempty(fnc)
                    fprintf('Index %d chosen\n', idx);
                else
                    feval(fnc, idx);
                end
            end
        end
        
        function CheckForUpdate(userIsAsking, app, appName)
            if nargin<2
                app=BasicMap;
            end
            app.setAppVersion([], UMAP.VERSION, SuhEpp.VERSION);
            try
                if nargin>2
                    temp=app.appName;
                    app.appName=appName;
                end
                SuhWebUpdate.Check(app, userIsAsking, .25, [], [], ...
                    ArgumentClinic.PROP_UPDATE_CHECK_TIME, 36);
                if nargin>2
                    app.appName=temp;
                end
            catch ex
                ex.reportProblem
            end
        end
        
        function jobWatch=RunJobs(folder)
            if ~isempty(folder)
                jobWatch=SuhJob(folder, [], false, ... %accept all args
                    @(this, job)go(job), 'south', ...
                    'ArgumentClinic.JobWatch', ...
                    {'suh_pipelines', 'run_umap', ...
                    'run_epp', 'run_match'});
            else
                jobWatch=[];
            end
            
            function finalArgs=go(job)
                finalArgs=[];
                if isfield(job, 'command')
                    cmd=job.command;
                    isPipeArgNeeded=true;
                    if strcmp(cmd, 'run_umap')
                        pipe='umap';
                    elseif strcmp(cmd, 'run_epp')
                        pipe='epp';
                    elseif strcmp(cmd, 'run_match')
                        pipe='match';
                    else %suh_pipelines
                        pipe=Args.GetIfNotPairs(...
                            'pipeline', job.args{:});
                        isPipeArgNeeded=false;
                        if isempty(pipe)
                            pipe='epp';%epp by default
                        end
                    end
                    if strcmpi(pipe, 'umap')...
                            || strcmp(pipe, 'epp')
                        job.props.set('data', job.args{1});
                        job.varArgs=job.args(2:end);
                    else
                        job.varArgs=job.args;
                    end
                    job.props.set('pipeline', pipe);  
                    if isPipeArgNeeded
                        job.varArgs=[job.varArgs 'pipe', pipe];
                    end
                end
                csv=job.props.get('data');
                if ~isempty(csv)
                    varArgs=Args.RemoveArg(job.varArgs, 'data');
                    obj=suh_pipelines(csv, varArgs{:});
                else
                    pipe=job.props.get('pipeline');
                    if isempty(pipe)
                        pipe=job.props.get('pipe');
                    end
                    if isempty(pipe)||~strcmpi(pipe,'match')
                        warn=sprintf(['<html>Received '...
                            'job with no data for pipeline '...
                            '%s in file %s<hr></html>'], ...
                            pipe, Html.FileTree(...
                            job.props.fileName));
                        msg(warn);
                        warning(['Received job with no data '...
                            'for pipeline %s in file %s'], pipe,...
                            job.props.fileName);                        
                        return;
                    else
                        varArgs=Args.RemoveArg(...
                            job.varArgs, 'pipeline');
                        varArgs=[...
                            'pipeline', pipe, varArgs];
                        obj=suh_pipelines(varArgs{:});
                    end
                end
                if ~isempty(obj)
                    finalArgs=obj.args;
                    finalArgs.done=true;
                end
            end
        end

        function fig=FlowJoBridgeFig(onlyIfOpen, reset)
            persistent singleton;
            if nargin>1 && reset
                if isempty(singleton)
                    delete(singleton);
                end
                singleton=[];
                fig=[];
                return;
            end
            if ~isempty(singleton) && ishandle(singleton)
                fig=singleton;
                figure(singleton);
                return;
            elseif nargin>0 && onlyIfOpen  
                fig=[];
                return;
            end
            singleton=ArgumentClinic.InitFig;
            fig=singleton;
        end

        function fig=InitFig
            app=BasicMap.Global;
            [fig, ~, personalized]=Gui.Figure(true, ...
                'FlowJoBride.Window', app, ...
                'south west', false);
            priorCloseFcn=get(fig, 'CloseRequestFcn');
            set(fig, 'Name', ['v' ArgumentClinic.VERSION ...
                ' Open bridge/pipelines'], ...
                'CloseRequestFcn', @(h,e)hush(h))
            %set(fig, 'visible', 'on')
            if ~personalized
                p=get(fig, 'position');
                p(3)=floor(p(3)/3.8);
                p(4)=floor(p(4)/4);
                app=BasicMap.Global;
                p(3)=app.adjustHighDef(p(3), 30);
                p(4)=app.adjustHighDef(p(4), 20);
                set(fig, 'position', p);
            end
            [~,bp]=Gui.ImageLabel(['<html>' ...
                '<center>Open v' num2str(ArgumentClinic.VERSION) ...
                ' of ....<hr>' ...
                Html.ImgXy('flowJoTM.png',[], .394, false,false, app) ...
                 '&nbsp;bridge<br></html>'], [], ArgumentClinic.FlowJoStory, ...
                @(h,e)flowJoBridge(h));
            dim=bp.getPreferredSize;
            dim.width=dim.width*1.5;
            dim.height=dim.height*1.5;
            bp.setPreferredSize(dim);
            h=Gui.SetJavaInFig(.02, .02, bp, ...
                Gui.LIGHT_YELLOW_COLOR, fig, app, false);
            set(h, 'units', 'normalized');
            set(h, 'Position', [0.001 0.1 1 1]);
            [~,bp]=Gui.ImageLabel( ...
                Html.WrapSmall('...or Herzenberg pipelines'), ...
                [], ArgumentClinic.FlowJoStory, ...
                @(h,e)ArgumentClinic());
            dim=bp.getPreferredSize;
            dim.width=dim.width*1.5;
            dim.height=dim.height*1.5;
            bp.setPreferredSize(dim);
            h=Gui.SetJavaInFig(.02, .02, bp, ...
                Gui.LIGHT_GRAY_COLOR, fig, app, false);
            set(h, 'units', 'normalized');
            set(h, 'Position', [0.001 0.001 1 1]);
            
            function flowJoBridge(h)
                FlowJoTree.Open(Gui.WindowAncestor(h));
            end

            function hush(h)
                try
                    feval(priorCloseFcn, h,[]);
                    [fjbCnt, fjbTrees]=FlowJoTrees;
                    if fjbCnt>0
                        fjbCnt=0;
                        N=length(fjbTrees);
                        for i=1:N
                            fjb=fjbTrees{i};
                            if isa(fjb, 'FlowJoTree') && ...
                                    Gui.IsVisible(fjb.fig)
                                fjbCnt=fjbCnt+1;
                            end
                        end
                    end
                    if fjbCnt>0
                        jw=Gui.JWindow(h);
                        str=String.Pluralize2('FlowJoBridge window', fjbCnt);
                        Gui.Wind;
                        askYesOrNo(struct( ...
                            'modal', false, ...
                            'checkFnc', @(jd, finalAnsw)ArgumentClinic.CheckClose(finalAnsw),...
                            'javaWindow', jw, ...
                            'msg', ['Close ' str '?']), 'Confirm...', ...
                            'north', false, '', 'FlowJoBridge.CloseAll');
                        delete(h);
                    else
                        delete(h);
                    end
                    app.save;
                catch ex
                    ex.getReport
                end
            end
        end

        function ok=CheckClose(finalAnsw)
            ok=true;
            if strcmpi('Yes', finalAnsw)
                try
                    close all
                catch ex
                    ex.getReport
                end
            end
        end
    end
end
