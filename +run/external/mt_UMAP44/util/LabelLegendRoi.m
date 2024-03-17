classdef LabelLegendRoi < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(SetAccess=private)
        xyData;
        javaLegend;
        btnCreate;
        btns;
        btnLbls;
        map;
        asked=false;
        buildingMultiplePolygons=false;
        dfltFile;
        save_output=false;
        output_folder=[];
        ax;
        roiCallback
        sprv;
        fncRefreshUst;
        fncRefreshMatch;
        btnUst;
        btnMlp;
        btnMatchOverlap;
        btnMatchSimilarity;
        pnlMatch;
        unreducedData;
        busy;
        jwMlpTip;
        args;
        useClusterBoundaries=true;
        suhRunUmap=[];
    end
    
    properties
        labels;    
        percentClosest=.94;
        props;
    end
    
    methods
        function this=LabelLegendRoi(xyData, roiCallback, ...
                javaLegend, btns, btnLbls, labels, map, args, suhRunUmap)
            if nargin>8
                this.suhRunUmap=suhRunUmap;
            end
            this.xyData=xyData;
            this.roiCallback=roiCallback;
            this.javaLegend=javaLegend;
            this.btns=btns;
            this.btnLbls=btnLbls;
            if isempty(map)
                this.map=Map;
            else
                this.map=map;
            end
            this.labels=labels;
            this.props=BasicMap.Global;
            if nargin>6
                this.args=args;
            end
        end
        
        function followParentFig(this, parentFig, where, close)
            SuhWindow.Follow(this.javaLegend, parentFig, where, close);
        end
        
        function [lbls, btns]=getSelected(this)
            [lbls, btns]=QfTable.GetSelectedLabels(this.btns, this.btnLbls);
        end
        
        
        function processSelections(this, grouped)
            assert(size(this.xyData, 1)==length(this.labels));
            [lbls, btns_]=this.getSelected;
            if grouped
                [roi, name, key]=this.makeOrDeleteOne(lbls);
                saveRoi;
                drawnow;
            else
                nChoices=length(lbls);
                fig=[];
                if ~isempty(this.suhRunUmap)
                    fig=this.suhRunUmap.fig;
                    if ~isempty(fig)
                        txt=String.Pluralize2('label choice', nChoices');
                        hBusy=Gui.ShowBusy(fig, Gui.YellowSmall('Creating polygon gates'));
                        ttl=sprintf('Converting %s to polygons', txt);
                        hWait=Gui.ProgressBar(ttl, true, fig, 'south++');
                        jBusy=Gui.ShowBusy(this.javaLegend, Gui.YellowSmall(ttl));
                    end
                end
                for i=1:nChoices
                    if ~isempty(fig)
                        stop=getappdata(hWait,'canceling');
                        if stop
                            break;
                        end
                    end
                    this.buildingMultiplePolygons=i<nChoices;
                    if ~this.useClusterBoundaries
                        [roi,name, key]=this.makeOrDeleteOne(lbls(i), btns_{i});
                        saveRoi;
                    else
                        this.createPolygonsForMajorityClustersOfLabel( ...
                            lbls(i), btns_{i});
                    end
                    if ~isempty(fig)
                        waitbar2a(i/nChoices, hWait, ...
                            sprintf('%d of %s done', i, txt));
                    end
                    drawnow;
                end
                if ~isempty(fig)
                    Gui.HideBusy(fig, hBusy);
                    Gui.HideBusy(this.javaLegend, jBusy);
                    delete(hWait);
                end
            end
            
            function saveRoi
                if ~isempty(roi) && ~isempty(this.args.save_roi)
                    try
                        if isempty(name)
                            name=LabelLegendRoi.Rename(this.map, roi);
                        end
                        if isempty(name)
                            delete(roi);
                            removed=this.map.remove(key);
                        else
                            feval(this.args.save_roi, key, ...
                                roi, name, this.xyData, this.args)
                        end
                    catch ex
                        ex.getReport
                    end
                end
            end
        end
        
        
        function [roi, name, key]=makeOrDeleteOne(this, lbls, btn)
            key=num2str(sort(lbls));
            name=[];
            if this.map.containsKey(key)
                if isempty(this.args.save_roi)
                    roi=this.map.remove(key);
                    if ~isempty(roi)
                        delete(roi);
                        this.map.remove([key '.name']);
                    end
                else
                    msg('Clusters already gated...', 0, 'south');
                    roi=[];
                end
                return;
            end
            if length(lbls)>1
                idxs=ismember(this.labels, lbls);
            else
                idxs=this.labels==lbls;
            end
            clump=this.xyData(idxs,:);
            [~,idxs]=MatBasics.PercentClosest(clump, this.percentClosest);
            [roi, ~]=RoiUtil.NewForXy(...
                this.ax, clump(idxs,:), .00121, this.roiCallback);
            this.map.set(key, roi);
            N=length(lbls);
            if N>1
                [name, cancelled]=inputDlg(...
                    'Enter summary name for properties', ...
                    'ROI name', key);
                if cancelled
                    name=key;
                end
            else
                name=LabelLegendRoi.Name(btn);
            end
            this.map.set([key '.name'], name);
        end
        
        function save(this)
            if LabelLegendRoi.Save(this.map, this.save_output,...
                    this.output_folder, this.dfltFile, ...
                    this.asked, reduction, this.args, this.suhRunUmap.sprv)
                this.asked=true;
            end
        end
        
        function createPolygonsForMajorityClustersOfLabel(this, label, btn)
            key=num2str(label);
            name=LabelLegendRoi.Name(btn);
            clues=this.suhRunUmap.clusterIdentifiers;
            u=unique(clues);
            u(u==0)=[];
            nClues=length(u);
            clustersForLabel=[];
            for i=1:nClues
                clue=u(i);
                labelsInCluster=this.labels(clues==clue);
                [cnts, descendingOrder, lbls]=MatBasics.HistCounts(labelsInCluster);
                if ~isempty(lbls)
                    top=descendingOrder(1);
                    if lbls(top)==label
                        fprintf(['Cluster %d contains %d instances ' ...
                            'of %s (%d)\n'], clue, cnts(top), name, label);
                        clustersForLabel(end+1)=clue;
                    end
                end
            end
            used=length(clustersForLabel);
            for i=1:used
                k=[key '.' num2str(i)];
                n=[name ' #' num2str(i)];
                if this.map.containsKey(k)
                    if isempty(this.args.save_roi)
                        roi=this.map.remove(k);
                        if ~isempty(roi)
                            delete(roi);
                            this.map.remove([k '.name']);
                        end
                    else
                        fprintf('Cluster %s already gated...\n', n);
                        
                    end
                    continue;
                end
                clue=clustersForLabel(i);
                clump=this.xyData(clues==clue,:);
                [roi, ~]=RoiUtil.NewForXy(...
                    this.ax, clump, .00121, this.roiCallback);
                this.map.set(k, roi);
                this.map.set([k '.name'], n);
                if ~isempty(this.args.save_roi)
                    feval(this.args.save_roi, k, ...
                        roi, n, this.xyData, this.args);
                end
            end
        end
        function updateLegendGui(this, ax, showCreateBtn, parentFig, where, close)
            this.ax=ax;
            if showCreateBtn && size(this.xyData, 2)==2
                jcmp=this.javaLegend.getContentPane;
                bp=Gui.BorderPanel([],2,2, 'Center', ...
                    Gui.Panel(Gui.NewBtn(...
                    Html.WrapSmallBold('Regions of interest'), ...
                    @(h,e)dropDownMenu(this, h), ...
                    'Create/remove polygons for above subset selection', ...
                    'polygonGate.png')));
                if ~isempty(this.btnMlp) ...
                        || ~isempty(this.pnlMatch)
                    west=Gui.FlowLeftPanel(1,0);
                    if ~isempty(this.btnMlp)
                        west.add(this.btnMlp)
                    end
                    if ~isempty(this.pnlMatch)
                        west.add(this.pnlMatch);
                    end
                    bp.add(west, 'West')
                end
                if ~isempty(this.btnUst)
                    bp.add(this.btnUst, 'East')
                end
                jcmp.add(bp, 'South')
                this.javaLegend.pack;
            end
            drawnow;
            this.followParentFig(parentFig, where, close)
        end
        
        function shakeMatch(this)
            if ~isempty(this.pnlMatch)
                edu.stanford.facs.swing.Basics.Shake(...
                    this.pnlMatch,5);
                MatBasics.RunLater(@(h,e)tip(), 2);
            end
            
            function tip
                BasicMap.Global.showToolTip(...
                    this.pnlMatch, ...
                    ['<html>You might ALSO wish to '...
                    '<br>refresh the match with<br>'...
                    'the prior classification??<hr>'...
                    '</html>'], -25, -40)

            end
        end
       
        function setSaveInfo(this, dfltFile, save_output, output_folder)
            this.dfltFile=dfltFile;
            this.save_output=save_output;
            this.output_folder=output_folder;
        end
        
        function resetFromPrior(this, prior)
            location=prior.javaLegend.getLocation;
            this.javaLegend.setLocation(location);
            Gui.SetJavaVisible(this.javaLegend);
            prior.javaLegend.dispose;
            this.shakeMatch;
            this.showMlpUstTip;
        end
    end
    
    methods(Access=private)
        
        function dropDownMenu(this, h)
            jm=PopUp.Menu;
            app_=BasicMap.Global;
            [lbls, btns_]=this.getSelected;
            nSelectedLbls=length(lbls);
            c={};
            enable=false;
            c{1}=Gui.NewMenuLabel(jm,'ROI options (regions of interest)');
            jm.addSeparator;
            callbackForIndividualOnes=...
                @(h,e)processSelections(this, false);
            if isempty(this.args.save_roi)
                verb='Create/remove';
            else
                verb='Create';
            end
            if nSelectedLbls>1
                c{2}=Gui.NewMenuItem(jm, sprintf(['%s individual'...
                    ' polygons for above %d selections'],...
                    verb, nSelectedLbls), callbackForIndividualOnes);
                enable=true;
            else
                if nSelectedLbls==0
                    word='NOTHING selected!!';
                else
                    word=Html.StripHtmlWord(char(btns_{1}.getText));
                end
                c{2}=Gui.NewMenuItem(jm,...
                    sprintf('<html>%s polygon(s) for "%s"</html>', ...
                    verb, word), callbackForIndividualOnes, [],[],...
                    nSelectedLbls==1);
            end
            c{end+1}=Gui.NewMenuItem(jm, sprintf(...
                '%s one polygon grouping all %d selections', ...
                verb, nSelectedLbls), ...
                @(h,e)processSelections(this, true),[],[],enable);
            c{end+1}=Gui.NewMenuItem(jm, ...
                'Save polygon information', @(h,e)saveLegendRois());
            if ~isempty(this.suhRunUmap.clusterIdentifiers)
                jm.addSeparator;
                c{end+1}=Gui.NewCheckBoxMenuItem(jm, ...
                    'Use cluster boundaries', @use, [], ...
                    Html.Wrap(['Polygons encompass each cluster that' ...
                    ' mostly consists <br>of a selection''s data points.  ' ...
                    'Otherwise the polygons <br>include all ' ...
                    'of a selection''s data points <b>regardless</b> of ' ...
                    '<br>cluster <u>boundaries</u>.']), this.useClusterBoundaries);
            end
            jm.show(h, 15, 15);
            
            function use(cb, e)
                this.useClusterBoundaries=cb.isSelected;
            end
           
            function saveLegendRois
                if this.map.size<1
                    msg('Nothing to save ...');
                else
                    this.save;
                end
            end
        end
    end

    methods

        function ok=setUst(this, sprv, fncRefreshUst, ...
                unreducedData, fncRefreshMatch)
            sprv.storeClusterIds=true;
            this.sprv=sprv;
            this.unreducedData=unreducedData;
            this.fncRefreshUst=fncRefreshUst;
            if nargin>4
                this.fncRefreshMatch=fncRefreshMatch;
            end
            img='ust.png';
            this.btnUst=Gui.ImageButton(img, ...
                'Manage classification via supervisors', @(h,e)showOptions(h));
            if ~isempty(sprv.mlp_confidence)
                this.btnMlp=Gui.ImageButton('mlp.png', ...
                    'Manage MLP confidence', @(h,e)showMlpUstTip(this));
            end
            if nargin>4 && ~isempty(fncRefreshMatch)
                this.btnMatchOverlap...
                    =Gui.ImageButton(...
                    'match16.png', ...
                    ['Match using the F-measure '...
                    'harmonic mean of overlap in/out'], ...
                    @(h,e)refreshMatch(this, 2));
                this.btnMatchSimilarity...
                    =Gui.ImageButton(...
                    'match.png', ...
                    ['Match using mass+distance'...
                    ' similarity (QFMatch)'], ...
                    @(h,e)refreshMatch(this, 1));
                this.pnlMatch=Gui.FlowLeftPanel(1, 0, ...
                    this.btnMatchOverlap, ...
                    this.btnMatchSimilarity);
            end

            function showOptions(h)
                jMenu=PopUp.Menu;
                this.addUstMenus(jMenu);
                jMenu.show(h, 5, 5);
            end
        end

        function refreshUst(this, pu)
            try
                if nargin<2
                    feval(this.fncRefreshUst, this.xyData, this.sprv.matchType, [], false);
                else
                    feval(this.fncRefreshUst, this.xyData, ...
                        this.sprv.matchType, pu, false);
                end
            catch ex
                ex.getReport
            end
        end

        function refreshMatch(this, matchStrategy)
            this.showBusy('Running Darya''s QFMatch', 'orlova.png', .33);
            feval(this.fncRefreshMatch, matchStrategy);
            this.hideBusy;
        end

        function matchType=getUstMatchType(this)
            matchType=this.props.getNumeric('ustMatchType', 3);
        end
        
        function pu=showBusyReclustering(this)
            this.showBusy('Re-running Guenther''s<br>clustering', 'guenther.png', .3);
        end

        function setUstMatchType(this, matchType)
            spr=this.sprv;
            [pu, txt]=Supervisors.StartUstRematch(...
                matchType, this.javaLegend);
            this.showBusy(txt, 'connor.png', .243);
            this.props.set('ustMatchType', num2str(matchType));
            spr.setMatchType(matchType);
            if matchType==4
                spr.computeNearestNeighborsUnreduced(...
                    this.unreducedData, pu);
            %elseif matchType>=2
            %    spr.computeNearestNeighbors(this.xyData, pu);
            end
            %if isempty(spr.lastClusterIds)
            %    [~, ~, spr.density]=spr.findClusters(this.xyData, pu);
            %end
            %spr.matchClusters(this.xyData(:,1:2), spr.density, ...
            %    spr.density.numClusters, spr.lastClusterIds, matchType, pu);
            this.refreshUst(pu);
            this.hideBusy;
            pu.close;
        end

        function addUstMenus(this, jMenu)
            spr=this.sprv;
            matchType=spr.matchType;
            mnu=jMenu;
            Gui.NewMenuLabel(mnu, 'UMAP Supervisor classification');
            mnu.addSeparator;
            %mnu=Gui.NewMenu(jMenu, 'Match UMAP template''s training sets with ...', ...
            %    [], 'match.png', false);
            
            ustM=UMAP.MATCH;
            ustTip=UMAP.MATCH_TIP;
            for j=1:UMAP.MATCH_N
                Gui.NewCheckBoxMenuItem(mnu, ustM{j},  ...
                    @(h,e)setUstMatchType(this, j-1), [],...
                    ustTip{j}, matchType==j-1);
            end
            mnu.addSeparator;
            Gui.NewMenuLabel(mnu, 'MLP classification');
            %this.addMlpMenuItems(mnu);
            this.addMlpMenuItems(jMenu, 'mlp.png');
            mnu.addSeparator;
            Gui.NewMenuLabel(mnu, 'Prior classification');
            enable=~isempty(this.fncRefreshMatch);
            Gui.NewMenuItem(mnu, 'Match by overlap',...
                @(h,e)refreshMatch(this, 2), 'match16.png',...
                'Match using the F-measure harmonic mean of overlap in/out', ...
                enable);
            Gui.NewMenuItem(mnu, 'Match by similarity', ...
                @(h,e)refreshMatch(this, 1), 'match.png',...
                'Match using mass+distance similarity', ...
                enable);
        end

        function addMlpMenuItems(this, mnu, icon)
            if nargin<3
                icon=[];
            end
            spr=this.sprv;
            if isempty(spr.mlp_use_python)
                mlpType='';
            elseif spr.mlp_use_python
                mlpType='MLP (TensorFlow)';
            else
                mlpType='MLP (fitcnet)';
            end
            if spr.mlp_ignore
                txt=['Use ' mlpType ' classifications'];
                tip=['Use ' mlpType ...
                    ' classification and UMAP overrides'];
            else
                txt=['Ignore ' mlpType ' classifications'];
                tip='Only use UMAP classification';
            end
            Gui.NewMenuItem(mnu, txt, ...
                @(h,e)toggleMlp(this), icon,...
                tip, ...
                ~isempty(spr.mlp_labels));
            
            hasConfidence=~isempty(spr.mlp_confidence) && ~spr.mlp_ignore;
            Gui.NewMenuItem(mnu, 'Ignore UMAP''s overrides of MLP',  ...
                @(h,e)setUstMlpConfidence(this, true),...
                'ust.png',...
                'Set MLP confidence level to 0 so UMAP overrides NOTHING', ...
                ~spr.mlp_ignore && hasConfidence && spr.mlp_confidence_level>0);
            Gui.NewMenuItem(mnu, sprintf('Explain UMAP''s %s MLP overrides', ...
                spr.mlp_overridden_txt),  ...
                @(h,e)showMlpUstTip(this), 'help2.png',...
                'See MLP confidence and oveerrides', ...
                hasConfidence);
            Gui.NewMenuItem(mnu, sprintf(...
                'Adjust override confidence level (%s)', ...
                String.encodePercent(...
                spr.mlp_confidence_level)),  ...
                @(h,e)setUstMlpConfidence(this),...
                'tool_data_cursor.gif',...
                'Adjust MLP confidence level', ...
                hasConfidence);
        end

        function toggleMlp(this)
            this.sprv.mlp_ignore=~this.sprv.mlp_ignore;
            preamble=['Classifying with UMAP''s  '...
                'supervisors<br>and '];
            if this.sprv.mlp_ignore
                this.showBusy([preamble ...
                    'UMAP methods only'],...
                    'connor.png', .21);
            else
                this.showBusy([preamble ' ', ...
                    'MLP+UMAP methods']);
            end
            this.refreshUst();
            this.hideBusy;
        end

        function level=setUstMlpConfidence(this, toZero)
            if nargin<2
                toZero=false;
            end
            spr=this.sprv;
            if toZero
                level=0;
                this.showBusy(['Ignoring UMAP overrides'...
                    ' with<br>MLP confidence level of 0%!!']);
            else
                this.showBusy;
                level=Mlp.GetConfidence(false, ...%for UST
                    [],  spr.mlp_confidence_level, ...
                    true, this.btnMlp, 'center');
            end
            if ~isempty(level)
                if level ~= spr.mlp_confidence_level
                    spr.setConfidenceLevel(level);
                    this.refreshUst();
                end
            end
            this.hideBusy;
        end
        
        function showMlpUstTip(this)
            if ~isempty(this.jwMlpTip) && this.jwMlpTip.isVisible
                this.jwMlpTip.requestFocus;
            else
                this.jwMlpTip=Supervisors.ShowMlpUstTip(this.sprv, ...
                    this.javaLegend, this.ax, this.xyData, ...
                    @reAdjust);
            end
            
            function reAdjust(h)
                newLvl=this.setUstMlpConfidence();
                if isempty(newLvl)
                    return;
                end
            end    
        end
       
        function showBusy(this, txt, icon, scale, showBusy)
            if nargin<5
                showBusy=true;
                if nargin<4
                    scale=.46;
                    if nargin<3
                        icon='ebrahimian.png';
                        if nargin<2
                            txt='Adjusting MLP confidence level';
                        end
                    end
                end
            end
            if ~isempty(this.jwMlpTip)
                this.jwMlpTip.dispose;
                this.jwMlpTip=[];
            end
            if BasicMap.Global.highDef
                scale=scale*1.9;
            end
            this.busy=Gui.ShowBusy(this.javaLegend, ...
                Gui.YellowSmall(txt, BasicMap.Global), icon, ...
                scale, showBusy);
        end
        
        function hideBusy(this)
            if ~isempty(this.busy)
                Gui.HideBusy(this.javaLegend, ...
                    this.busy, true);
                this.busy=[];
            end
        end

    end
    
    methods(Static)
        function name=Name(btn)
            app=BasicMap.Global;
            name=char(btn.getText);
            idx=String.IndexOf(name, app.supStart);
            if idx>0
                name=name(1:idx-1);
            end
            name=Html.Remove(name);
            name=strrep(name, '&bull;', '');
        end
        
        function ok=Save(map, save_output, output_folder, dfltFile, ...
                asked, reduction, args, oneNameOnly, sprv)
            if nargin<9
                sprv=[];
                if nargin<8
                    oneNameOnly='';
                end
            end
            ok=false;
            if map.size<1
                msg('No named regions to save yet...');
                return;
            end
            if isempty(args.save_roi)
                if ~save_output
                    [folder, file]=uiPutFile(output_folder, dfltFile, ...
                        BasicMap.Global, 'LabelLegendRoi.save', ...
                        'Save labeled ROI properties');
                    if isempty(folder)
                        return;
                    end
                    file=fullfile(folder,file);
                else
                    file=fullfile(output_folder, dfltFile);
                end
                props=java.util.Properties;
                keys=map.keys;
                N=length(keys);
                for i=1:N
                    key=keys{i};
                    if ~endsWith(key, '.name') && ~strcmp('keys', key)
                        roi=map.get(key);
                        str=MatBasics.XyToString(...
                            RoiUtil.Position(roi));
                        name=map.get([key '.name']);
                        props.setProperty([key '.position'],  str);
                        props.setProperty([key '.name'],  name);
                        props.setProperty([key '.type'],  class(roi));
                    end
                end
                fldr=fileparts(file);
                File.mkDir(fldr);
                File.SaveProperties2(file, props);
                if ~asked
                    File.OpenFolderWindow(file, '', true);
                end
            else
                keys=map.keys;
                N=length(keys);
                cnt=0;
                for i=1:N
                    key=keys{i};
                    if ~endsWith(key, '.name') && ~strcmp('keys', key)
                        cnt=cnt+1;
                    end
                end
                done=0;
                for i=1:N
                    key=keys{i};
                    if ~endsWith(key, '.name') && ~strcmp('keys', key)
                        roi=map.get(key);
                        name=map.get([key '.name']);
                        try
                            done=done+1;
                            proceed=isempty(oneNameOnly) ...
                                || strcmp(name, oneNameOnly);
                            if proceed
                                try
                                    feval(args.save_roi, key, roi, name, ...
                                        reduction, args, true, sprv)
                                catch
                                    feval(args.save_roi, key, roi, name, ...
                                        reduction, args)
                                end
                            end
                        catch ex
                            Gui.MsgException(ex, 'Gate saving issue?', 11);
                        end
                    end
                end
            end
            ok=true;
        end
        
        function Delete(map, roi)
            [key, name]=LabelLegendRoi.Find(map,roi);
            if ~isempty(key)
                roi=map.remove(key);
                if ~isempty(roi)
                    delete(roi);
                    map.remove([key '.name']);
                end
                fprintf('Deleted ROI, key=%s, name=%s]n', key, name);
            else
                % not identified with key and name yet
                delete(roi);
            end
        end
        
        function name=Rename(map, roi, suggestedName, tip)
            name='';
            if ~isvalid(roi)
                return;
            end
            [key,name]=LabelLegendRoi.Find(map,roi);
            if isempty(key)
                isNew=true;
                key=map.get('keys');
                if isempty(key)
                    key=1;
                end
                map.set('keys', key+1);
                key=num2str(0-key);
                map.set(key, roi);
            else
                isNew=false;
            end
            if isempty(name)
                if nargin>2
                    name=suggestedName;
                end
            end
            oldClr=RoiUtil.SetColor(roi, 'red');
            if nargin>3
                ttl=tip;
            else
                ttl='';
            end
            while true
                [name, cancelled]=inputDlg( ...
                    struct('where', 'north+', 'choiceTitle', ttl, ...
                    'msg', 'Name this region of interest (ROI)'), ...
                    'ROI names', name);
                if cancelled
                    break;
                end
                otherKey=LabelLegendRoi.FindName(map, name);
                if ~isempty(otherKey) && ~isequal(otherKey, key)
                    msgWarning(Html.SprintfHr(['The name "<b>%s</b>" ' ...
                        '<br>ALREADY exists ... enter a new name!'], name), 5)
                else
                    break;
                end
            end
            RoiUtil.SetColor(roi, oldClr);
            if ~cancelled
                RoiUtil.SetLabel(roi, name);
                map.set([key '.name'], name);
                if isNew
                    map.set(key, roi);
                end
            end
        end
        
        function [key, name]=Find(map, roi)
            keys=map.map.keys;
            rois=map.map.values;
            N=length(rois);
            for i=1:N
                if isequal(roi, rois{i})
                    key=keys{i};
                    name=map.get([key '.name']);
                    return;
                end
            end
            key='';
            name='';
        end
        
        function key=FindName(map, name)
            keys=map.map.keys;
            N=length(keys);
            for i=1:N
                key=keys{i};
                name2=map.get([key '.name']);
                if strcmpi(name, name2)
                    return;
                end
            end
            key='';
        end
        
    end
end