classdef Plots < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(Constant)
        DOING_ALL='doingAll';
        DEBUGGING=false;
    end
    
    properties(SetAccess=private)
        is3D;
        legendHs;
        mdns;
        stds;
        legendH;
        otherPlots={};
    end
    
    properties
        Hs;
        otherHs;
        cnts;
        names;
        N;
        distFactor;
        mx;
        mxI;
        mn;
        mnI;
        ax;
        javaLegend;
        otherPlotMap;
    end
    
    methods
        function this=Plots()
            this.Hs=[];
        end
        
        function setHs(this, Hs)
            this.Hs=Hs;
            this.N=length(this.Hs);
            this.cnts=zeros(1, this.N);
            if isempty(this.names) 
                this.names=cell(1, this.N);
                for i=1:this.N
                    this.cnts(i)=length(get(this.Hs(i), 'XData'));
                    this.names{i}=get(this.Hs(i), 'DisplayName');
                end
            else
                for i=1:this.N
                    this.cnts(i)=length(get(this.Hs(i), 'XData'));
                end
            end
            this.setMinMax;
            this.ax=get(Hs(1), 'Parent');
        end
        
        function setMinMax(this)
            [this.mx, this.mxI]=max(this.cnts);
            [this.mn, this.mnI]=min(this.cnts);
        end
        
        function setOtherHs(this, Hs)
            this.otherHs=Hs;
        end
        
        function setLegendHs(this, Hs)
            this.legendHs=Hs;
        end
        
        function initStats(this)
            if isempty(this.mdns)
                if isempty(this.Hs)
                    return;
                end
                [this.mdns, this.stds, this.is3D]=Plots.Stats(this.Hs);
                if this.is3D
                    this.distFactor=6;
                else
                    this.distFactor=3;
                end
            end
        end        
        
        function setNames(this, names)
            this.names=names;
        end
        
        function names=getNames(this)
            names=this.names;
        end
        
        function toggleVisibility(this, idx)
            if ~isempty(this.otherHs)
                R=size(this.otherHs,1);
                for r=1:R
                    plotH=this.otherHs(r, idx);
                    if strcmp('off', get(plotH, 'visible'))
                        set(plotH(r), 'visible', 'on');
                    else
                        set(plotH(r), 'visible', 'off');
                    end
                end
            end
        end
        
        function clear(this)
            try
                N_=length(this.Hs);
                for i=1:N_
                    delete(this.Hs(i));
                end
                N_=length(this.otherHs);
                for i=1:N_
                    delete(this.otherHs(i));
                end
            catch ex
                %ex.getReport
            end
        end
        
        function addOtherPlots(this, other)
            this.otherPlots{end+1}=other;
        end
        
        function flashOtherPlots(this, idx, on)
           if ~isempty(this.otherPlots)
               name=edu.stanford.facs.swing.Basics.RemoveXml(String.RemoveTex(this.names{idx}));
               lIdx=String.LastIndexOf(name, ' training ');
               if lIdx>0
                   name=name.substring(0,lIdx-1);
               end
               nOthers=length(this.otherPlots);
               for i=1:nOthers
                   other=this.otherPlots{i};
                   nNames=length(other.names);
                   for j=1:nNames
                       if strcmp(name, edu.stanford.facs.swing.Basics.RemoveXml(String.RemoveTex(other.names{j})))
                           perc=other.cnts(j)/sum(other.cnts);
                           otherH=other.Hs(j);
                           set(otherH, 'MarkerEdgeColor', ...
                               get(this.Hs(idx), 'MarkerEdgeColor'));
                           if nargin<3
                               if isequal(get(otherH, 'visible'), 'on')
                                   set(otherH, 'visible', 'off');
                                   times=3;
                               else
                                   set(otherH, 'visible', 'on');
                                   times=4;
                               end
                           else
                               if ~on
                                   set(otherH, 'visible', 'off');
                                   times=3;
                               else
                                   set(otherH, 'visible', 'on');
                                   times=4;
                               end
                           end
                           Plots.Flash(otherH, perc, times);
                           break;
                       end
                   end
               end
           end
           if nargin>2
               this.setOtherPlotMapVisible(idx, on)
           end
        end
        
        function yes=isInOtherPlotMap(this, idx)
            if ~isempty(this.otherPlotMap)
                H=this.Hs(idx);
                yes=this.otherPlotMap.isKey(H);
            else
                yes=false;
            end
        end
        
        function setOtherPlotMapVisible(this, idx, on)
            if ~isempty(this.otherPlotMap)
                H=this.Hs(idx);
                if on
                    vis='on';
                else
                    vis='off';
                end
                if this.otherPlotMap.isKey(H)
                    v=this.otherPlotMap(H);
                    if iscell(v)
                        v=v{1};
                    end
                    if any(v<=0)
                        v(v<0)=v;
                    end
                    set(v, 'visible', vis);
                end
            end
        end
        
        function setOtherPlotVisible(this, idx, on)
           if ~isempty(this.otherPlots)
               name=edu.stanford.facs.swing.Basics.RemoveXml(String.RemoveTex(this.names{idx}));
               lIdx=String.LastIndexOf(name, ' training ');
               if lIdx>0
                   name=name.substring(0,lIdx-1);
               end
               nOthers=length(this.otherPlots);
               for i=1:nOthers
                   other=this.otherPlots{i};
                   nNames=length(other.names);
                   for j=1:nNames
                       if strcmp(name, edu.stanford.facs.swing.Basics.RemoveXml(String.RemoveTex(other.names{j})))
                           perc=other.cnts(j)/sum(other.cnts);
                           otherH=other.Hs(j);
                           if ~on
                               set(otherH, 'visible', 'off'); 
                           else
                               set(otherH, 'visible', 'on');
                           end
                       end
                   end
               end
           end
           this.setOtherPlotMapVisible(idx, on)
        end
    end
    
    methods(Static)
        
        function [mdns, stds, is3D]=Stats(plotHs)
            N=length(plotHs);
            mdns=zeros(N,3);
            stds=zeros(N,3);
            is3D=true;
            for ii=1:N
                plotH=plotHs(ii);
                if ~ishandle(plotH)
                    return;
                end
                x_=get(plotH, 'XData');
                y_=get(plotH, 'YData') ;
                z_=get(plotH, 'ZData');
                if isempty(z_)
                    [RR,CC]=size(get(plotHs(ii), 'XData'));
                    z_=ones(RR, CC);
                    is3D=false;
                end
                d=[x_; y_; z_]';
                mdns(ii,:)=median(d);
                stds(ii,:)=std(d);
            end
        end
        
        function Flash(plotH, perc, times)
            vi=get(plotH, 'Visible');
            ms=get(plotH, 'MarkerSize');
            if perc<.02
                set(plotH, 'MarkerSize', ms+7)
            elseif perc<.05
                set(plotH, 'MarkerSize', ms+5)
            elseif perc<.1
                set(plotH, 'MarkerSize', ms+3);
            elseif perc<.2
                set(plotH, 'MarkerSize', ms+2);
            end
            for i=1:times
                %MatBasics.RunLater(@(h,e)go, .1);
                go
            end
            %MatBasics.RunLater(@(h,e)done, .51);
            done
            function go
                set(plotH, 'visible', 'off');
                pause(0.1);
                set(plotH, 'visible', 'on'); 
                pause(0.05);
            end
            function done
                set(plotH, 'visible', vi);
                set(plotH, 'MarkerSize', ms);
            end
        end
        
        
        function [HL, fcnMotion, btns, sortI, plots, sortGui]=Legend(...
                plotHs, names, idxsWithName, xMargin, yMargin, doMotion, ...
                freq, priorFcn, doubleClicker, doJavaLegend, oldJavaBtns,...
                southComponent, selectedCallback)
            sortGui=[];
            if nargin<13
                selectedCallback=[];
                if nargin<12
                    southComponent=[];
                    if nargin<11
                        oldJavaBtns=[];
                        if nargin<10
                            if nargin<9
                                doubleClicker=[];
                                if nargin<8
                                    priorFcn=[];
                                    if nargin<7
                                        freq=[];
                                        if nargin<6
                                            doMotion=true;
                                            if nargin<5
                                                yMargin=[];
                                                if nargin<4
                                                    xMargin=[];
                                                    if nargin<3
                                                        idxsWithName=[];
                                                        if nargin<2
                                                            names={};
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            doJavaLegend=false;
                        end
                    end
                end
            end
            app=BasicMap.Global;
            if doJavaLegend
                sup1=app.supStart;
                sup2=app.supEnd;
            end
            btns=[];
            if isa(plotHs, 'Plots')
                that=plotHs;
                plots=plotHs;
                N=plots.N;
            else
                that=[];
                plots=Plots;
                plots.setNames(names);
                if iscell(plotHs)
                    N=length(plotHs);
                    plotHs2=zeros(1,N);
                    for i=1:N
                        plotHs2(i)=plotHs{i};
                    end
                    plotHs=plotHs2;
                else
                    [R,N]=size(plotHs);
                    if R>1
                        otherPlotHs=plotHs(2:end,:);
                        plotHs=plotHs(1,:);
                        plots.setOtherHs(otherPlotHs);
                    end
                end
                plots.setHs(plotHs);
                names=plots.names;
            end
            if isempty(doJavaLegend)
                HL=[];
                fcnMotion=[];
                sortI=[];
                return;
            end
            if N<1
                return;
            end              
            if isempty(idxsWithName)
                idxsWithName=1:N;
            end
            assert(N>=length(idxsWithName), ...
                '# of plotHs must be >= # of idxsWithName');
            ax=get(plots.Hs(1), 'Parent');
            fig=get(ax, 'Parent');
            curOver=0;
            overH=[];
            try
                set(fig,'WindowButtonMotionFcn', @legMotion);
            catch
            end
            fcnMotion=@legMotion;
            X=xlim(ax);
            Y=ylim(ax);
            assert( length(names) == length(idxsWithName), ...
                'Must be as many namedPlotIdxs as names');

            mx=plots.mx;
            cnts=plots.cnts;
            if isempty(freq)
                Cnt=sum(cnts);
            else
                if length(freq)==length(cnts)
                    if any(isnan(freq))
                        freq(isnan(freq))=-1;
                    end
                    Cnt=sum(freq(freq>=0));
                    cnts=freq;
                    mx=max(cnts);
                else
                    Cnt=freq;
                end
            end
            nNames=length(names);
            cnts2=zeros(1, nNames);
            cntPref=app.getNumeric('cntPref', 2);
            if ~doJavaLegend
                legendHs=zeros(1, nNames);
                for i=1:nNames
                    idx=idxsWithName(i);
                    cnt=cnts(idx);
                    cnts2(i)=cnt;
                    if cnt>=0
                        strCnt=['    ^{\color{blue}' ...
                            String.encodePercent(cnt, Cnt, 1) '  '...
                            '\color[rgb]{.2 .5 .2} ' ...
                            String.encodeCount(cnt, cntPref) '} '];
                        mkSz=6 + (cnt/mx*15);
                        clr=get(plots.Hs(idx), 'MarkerEdgeColor');
                    
                        legendHs(i)=plot(ax, X(1)-X(2), Y(1)-Y(2), 's', ...
                            'MarkerSize', mkSz, 'Color',  clr, ...
                            'MarkerFaceColor', clr, 'LineStyle', 'none');
                    else
                        if isa(plotHs, 'Plots')
                            clr=get(plots.Hs(idx), 'MarkerEdgeColor');
                            legendHs(i)=plot(ax, X(1)-X(2), Y(1)-Y(2), 's', ...
                                'MarkerSize', 6, 'Color',  clr, ...
                                'MarkerFaceColor', clr, 'LineStyle', 'none');
                        else
                            legendHs(i)=plotHs(i);
                        end
                        names{i}=names{i};
                    end
                    
                end
                xlim(ax,X);
                ylim(ax,Y);
                [~,sortI]=sort(cnts2, 'descend');
                HL=legend(ax, legendHs(sortI), names(sortI), ...
                    'Location', 'northeast','AutoUpdate', 'off');
                if ~isempty(xMargin) && ~isempty(yMargin)
                    p=get(HL, 'position');
                    p2=Gui.GetNormalized(ax);
                    xNudge=1-(p2(1)+p2(3));
                    xNudge=xNudge-xMargin;
                    yNudge=1-(p2(2)+p2(4));
                    yNudge=yNudge-yMargin;
                    set(HL, 'position', [p(1)+xNudge, p(2)+yNudge p(3) p(4)]);
                end
                plots.setLegendHs(legendHs);
                HL.ItemHitFcn=@(h,event)Plots.LegendClick(plots, event,...
                    idxsWithName, cnts2);
            else
                for i=1:nNames
                    idx=idxsWithName(i);
                    cnt=cnts(idx);
                    cnts2(i)=cnt;
                    clr=get(plots.Hs(idx), 'MarkerEdgeColor');
                    if cnt>=0
                        strSymbol=['<font  ' Gui.HtmlHexColor(clr)...
                            '>&bull;</font>'];
                        f=11+(cnt/mx*15);
                        f= ceil((f-4)/4)+3;
                        if f>=10
                            f=10;
                        end
                        f=String.encodeInteger(f);
                        str=['<font size="' f '">' strSymbol '</font>'];
                        strCnt=[sup1 ' <font color="blue">' ...
                            String.encodePercent(cnt,Cnt,1) '</font>'...
                            '  <font ' Html.Color([.2 .5 .2]) ' ><i>' ...
                            String.encodeCount(cnt, cntPref) '</i></font>' ...
                            sup2 ];
                    else
                        ls=get(plots.Hs(idx), 'LineStyle');
                        if strcmpi(ls, 'none')
                            str=['<font size="7" ' Gui.HtmlHexColor(clr)...
                                '>&bull;</font>&nbsp;'];
                        else
                            str=['<font size="8" ' Gui.HtmlHexColor(clr)...
                                '>&minus;</font>&nbsp;'];
                        end
                        strCnt='';
                    end
                    names{i}=['<html>' str names{i} '&nbsp;&nbsp;'...
                        strCnt Html.EncodeSort('name', strtrim(char(...
                        edu.stanford.facs.swing.Basics.RemoveXml(...
                        lower(names{i}))))) ...
                        Html.EncodeSort('frequency', cnt) ...
                        '</html>'];
                end
                [~, sortI]=sort(cnts2, 'descend');
                [outerPnl, ~, btns,sortGui, innerPnl]=...
                    Radio.Panel3(names(sortI), 1:nNames, 13, ...
                    @(h,e, i)Plots.JavaLegendClick(plots, i, sortI,...
                    idxsWithName, plots.Hs, cnts2, h, selectedCallback), ...
                    true,Html.Wrap(['<i>Deselecting '...
                    '<b>hides</b> items in plot & selecting '...
                    '<b>unhides</b></i>']));
                innerPnl.setBackground(java.awt.Color(1, 1, 1))
                
                sortGui.setProperties(app, 'javaLegend.');
                if ~isempty(oldJavaBtns)
                    it=oldJavaBtns.iterator;
                    while it.hasNext
                        btn=it.next;
                        if ~btn.isSelected
                            k=StringArray.IndexOf(names,char(btn.getText));
                            m=find(sortI==k, 1);
                            if ~isempty(m)
                                b=btns.get(m-1);
                                b.setSelected(false);
                                set(plots.Hs(k), 'visible', 'off');
                                disp('was OFF');
                            end
                        end
                    end
                    sortGui.setAllChbText;
                elseif isequal(idxsWithName, 1:nNames)
                    for i=1:nNames
                        if isequal('off', get(plots.Hs(i), 'visible'))
                            b=btns.get(i-1);
                            b.setSelected(false);
                        end
                    end
                end
                if ~isempty(southComponent)
                    bp=Gui.BorderPanel(0,0);
                    bp.add(outerPnl, 'Center');
                    bp.add(southComponent, 'South');
                    outerPnl=bp;
                end
                if ~isempty(oldJavaBtns) && oldJavaBtns.size>0
                    HL=Gui.WindowAncestor(oldJavaBtns.get(0));
                    psz=HL.getPreferredSize;
                    HL.setContentPane(outerPnl);
                    HL.pack;
                    HL.setPreferredSize(psz);
                else
                    pu=showMsg(outerPnl, 'Legend', 'north east+', false, ...
                        false, 0, false);
                    HL=pu.dlg;
                end
                sortGui.dlg=HL;
                fncClose=get(fig, 'CloseRequestFcn');
                set(fig, 'CloseRequestFcn', @disposeJavaLegend);
                btns.get(0).getParent.scrollRectToVisible(java.awt.Rectangle(0,0,1,1))
                %if ismac && HL.isVisible && isempty(oldJavaBtns)
                %    if 1<size(get(0, 'MonitorPositions'), 1)
                %        if isequal('on', get(fig, 'visible'))
                %            MatBasics.RunLater(@(h,e)relocate(getjframe(fig)), .52)
                %        else
                %            disp('huh');
                %        end
                %    end
                %end
            end
            if ~isempty(that)
                that.legendH=HL;
            end
            
            function disposeJavaLegend(h,e)
                HL.dispose;
                if ischar(fncClose)
                    feval(fncClose);
                else
                    feval(fncClose, h, e)
                end
            end
            
            function relocate(ref)
                Gui.LocateJava(HL, ref, 'north east+');
                HL.setVisible(true);
            end
            function legMotion(hObj, event)
                if ~isvalid(ax)
                    return;
                end
                cp=get(ax, 'CurrentPoint');
                if isempty(cp)
                    return;
                end
                x=cp(2,1);
                y=cp(2,2);
                z=cp(2,3);
                try
                    e=Gui.GetPixels(HL);
                catch
                    e=[];
                end
                %[normX normY normFigX normFigY e]
                if ~isempty(e)
                    cp2=get(get(HL, 'Parent'), 'currentpoint');
                    if cp2(2)<=e(2)+e(4) && cp2(2)>=e(2)
                        if cp2(1)>=e(1) && cp2(1)<= e(1)+e(3)
                            if ~isempty(doubleClicker)
                                doubleClicker.stopListening;
                            end
                            return;
                        end
                    end
                end
                if ~isempty(doubleClicker)
                    doubleClicker.startListening;
                end
                if ~doMotion
                    if ~isempty(priorFcn)
                        feval(priorFcn, hObj, event);
                    end
                    return;
                end
                plots.initStats;
                if plots.is3D
                    [D,II]=pdist2(plots.mdns, [x y z], 'euclidean', 'Smallest', 1);
                    limit=pdist2(plots.mdns(II,:), ...
                        plots.mdns(II,:)+(plots.stds(II,:)*plots.distFactor));
                else
                    [D,II]=pdist2(plots.mdns(:,1:2), [x y], 'euclidean', 'Smallest', 1);
                    limit=pdist2(plots.mdns(II,1:2), plots.mdns(II,1:2)+...
                        (plots.stds(II,1:2)*plots.distFactor));
                end
                if Plots.DEBUGGING
                    [x y z ; II D limit] %#ok<NOPRT> 
                    plots.names{II} %#ok<NOPRT> 
                    if II==plots.mxI
                        disp('hey');
                    end
                end
                if D<=limit
                    over=II;
                else
                    over=0;
                end
                if ~ishandle(overH)
                    return;
                end
                if curOver ~= over
                    curOver=over;
                    if over>0                       
                        nameIdx=find(idxsWithName==II,1);
                        if nameIdx>0
                            nm=names{nameIdx};
                            if doJavaLegend
                                nm=strrep(nm, '<sup>', '_{{');
                                nm=strrep(nm, '</sup>', '}');
                                nm=char(...
                                    edu.stanford.facs.swing.Basics.RemoveXml(nm));
                                nm=strrep(nm, '_{{', '^{');
                                nm=strrep(nm, '&bull;','');                                
                            end
                            if isempty(overH)
                                if plots.is3D
                                    overH=text(ax, x, y, z, nm, ...
                                        'fontSize', 9, ...
                                        'color', [0 0 .5],...
                                        'EdgeColor', 'red', ...
                                        'FontName', 'Arial', ...
                                        'backgroundColor', [255/256 252/256 170/256]);
                                else
                                    overH=text(ax, x, y, nm, ...
                                        'fontSize', 9, 'color', [0 0 .5], 'EdgeColor',...
                                        'red', 'FontName', 'Arial', ...
                                        'backgroundColor', [255/256 252/256 170/256]);
                                end
                            else
                                if plots.is3D
                                    set(overH, 'visible', 'on', 'Position', ...
                                        [x y z], 'String', nm);
                                    %uistack(overH, 'top');
                                else
                                    set(overH, 'visible', 'on', 'Position', ...
                                        [x y 0], 'String', nm);
                                end
                            end
                        else
                            disp('background');
                        end
                    else
                        %disp('Over NOTHING');
                        if ~isempty(overH)
                            set(overH, 'visible', 'off');
                        end
                    end
                end
                if ~isempty(priorFcn)
                    feval(priorFcn, hObj, event);
                end
            end
        end
       
        function LegendClick(this, event, idxsWithName, nums)
            H=get(event.Peer, 'UserData');
            idx_=find(this.legendHs==event.Peer, 1);
            idx=idxsWithName(idx_);
            if ~isempty(H)
                set(event.Peer, 'UserData', []);
                delete(H);
                this.toggleVisibility(idx);
                return;
            end
            this.initStats;
            plotH=this.Hs(idx);
            perc= nums(idx_)/max(nums);
            Plots.Flash(plotH, perc, 4);
            X=this.mdns(idx, 1);
            Y=this.mdns(idx, 2);
            ax=get(this.Hs(idx), 'Parent');
            str=get(event.Peer, 'DisplayName');
            fsz=11;
            if this.is3D
                Z=this.mdns(idx, 3);
                Plots.ShowLegendTip(ax, event.Peer, X, Y, str, fsz, 0, Z);
                [X Y Z]
            else
                Plots.ShowLegendTip(ax, event.Peer, X, Y, str, fsz, 0);
            end
            this.toggleVisibility(idx);   
            this.flashOtherPlots(idx);            
        end
        
        function ok=JavaLegendClick(this, idx_, sortI, idxsWithName, plotHs, ...
                nums, h, selectedCallback)
            idx=idxsWithName(sortI(idx_));
            plotH=plotHs(idx);
            if ~ishandle(plotH)
                warning('Supervisors not showing');
                ok=false;
                return;
            end
            ok=true;
            hasOtherMap=this.isInOtherPlotMap(idx);
            doingAll=nargin==8 && isequal(Plots.DOING_ALL, ...
                    char(h.getActionCommand));
            if ~h.isSelected
                set(plotH, 'visible', 'off');
                times=3;
            else
                if doingAll && hasOtherMap
                    %nope
                    if this.otherPlotMap.size<length(this.Hs)
                        h.setSelected(false)
                        return;
                    end
                end
                set(plotH, 'visible', 'on');
                times=4;
            end
            if doingAll
                if ~isempty(selectedCallback)
                    feval(selectedCallback, h, idx, true);
                end
                this.setOtherPlotVisible(idx, h.isSelected);
                return;
            end
            perc=nums(idx)/max(nums);
            Plots.Flash(plotH, perc, times);
            if ~isempty(selectedCallback)
                feval(selectedCallback, h, idx, false);
            end
            this.flashOtherPlots(idx, h.isSelected);            
        end
        
        function ShowLegendTip(ax, hLegend, posX, posY, str, fsz, secs, posZ)
            if nargin<8
                [normX, normY]=Gui.DataToAxNorm(ax, posX, posY);
                tipH=text(normX+.03, normY+.03, str, ...
                    'fontSize', fsz, 'color', [0 0 .5], 'EdgeColor',...
                    'red', 'FontName', 'Arial', 'parent', ax, ...
                    'FontAngle', 'italic', 'units', 'normalized',...
                    'backgroundColor', [255/256 252/256 170/256]);
            else
                tipH=text(ax, posX, posY, posZ, str, ...
                    'fontSize', fsz, 'color', [0 0 .5], 'EdgeColor',...
                    'red', 'FontName', 'Arial', 'FontAngle', 'italic',...
                    'backgroundColor', [255/256 252/256 170/256]);
                
            end
            set(tipH, 'ButtonDownFcn', @(h,e)freeze(h, hLegend));
            if secs>0
                MatBasics.RunLater(@(h,e)closeTip, secs);
            else
                freeze(tipH, hLegend);
            end
            
            function freeze(hText, hLegend)
                if isempty(get(tipH, 'UserData'))
                    set(hText, 'UserData', true);
                    set(hText, 'FontSize', fsz-2);
                    set(hText, 'FontAngle', 'normal');
                    set(hLegend, 'UserData', hText);
                else
                    set(hLegend, 'UserData', []);
                    delete(hText);
                end
            end
            
            function closeTip
                if isempty(get(tipH, 'UserData'))
                    delete(tipH);
                end
            end
        end

        function fig=FreqScatterPlot(ConfusionMat, labels, figTtl, clrs)
            if nargin<4
                clrs=[];
            end
             N=length(labels);
             if isempty(clrs)
                 clrs=zeros(N,3);
                 for i=1:N
                     clrs(i,:)=Gui.HslColor(i,N);
                 end
             end
            if isnumeric(labels)
                lol=labels~=0;
            else
                lol=true(1, length(labels));
            end
            
            True_Freq = sum(ConfusionMat,2)./sum(sum(ConfusionMat));
            Predicted_Freq = sum(ConfusionMat,1)'./sum(sum(ConfusionMat));
            Max_Freq_diff = max(abs(True_Freq-Predicted_Freq))*100;
            disp(['delta_f = ' num2str(Max_Freq_diff)])
            %% Population Frequency scatter plot
            X=log(True_Freq*100);
            Y=log(Predicted_Freq*100);
            fig=Gui.NewFigure(true, false, true);
            ax=Gui.Axes(fig);
            scatter(ax, X(lol), Y(lol), 100, clrs(lol,:), 'filled')
            box(ax, 'on');
            grid(ax, 'on');
            xlabel(ax, 'Log(True frequency %)'),ylabel(ax, 'Log(Predicted frequency %)')
            title(ax, figTtl);
            xl=xlim(ax);
            yl=ylim(ax);
            [~,I]=sort(Y);
            flip=false(1, N);
            for i=1:N
                if mod(i,2)==0
                    flip(I(i))=true;
                end
            end
            xNudge=(xl(2)-xl(1))/50;
            yNudge=(yl(2)-yl(1))/50;
            for k=1:N
                if isnumeric(labels)
                    if labels(k)==0
                        continue
                    end
                    H=text(ax, X(k)+(xNudge),Y(k), num2str(labels(k)), ...
                        'FontSize', 7, 'FontWeight', 'normal');
                else
                    H=text(ax, X(k)+(xNudge),Y(k),labels{k}, ...
                        'FontSize', 7, 'FontWeight', 'normal');
                end
                if flip(k)
                    e=get(H, 'Extent');
                    p=get(H, 'Position');
                    p(1)=p(1)-(e(3)+(xNudge*2));
                    p(2)=p(2)+yNudge;
                    set(H, 'Position', p);
                end
            end
            lsline(ax);
            set(fig, 'Visible', 'on');
        end

        function fig=FreqBarPlot(ConfusionMat, labels, figTtl)
            if isnumeric(labels)
                lol=labels~=0;
            else
                lol=true(1, length(labels));
            end            
            True_Freq = sum(ConfusionMat,2)./sum(sum(ConfusionMat));
            Predicted_Freq = sum(ConfusionMat,1)'./sum(sum(ConfusionMat));
            Max_Freq_diff = max(abs(True_Freq-Predicted_Freq))*100;
            disp(['delta_f = ' num2str(Max_Freq_diff)])
            fig=Gui.NewFigure(true, false, true);
            ax=Gui.Axes(fig);
            [~,I]=sort(upper(labels));
            I=I(lol);
            bar(ax, [True_Freq(I)*100 Predicted_Freq(I)*100])
            xticks(ax, 1:24)
            xticklabels(ax, labels(I))
            xtickangle(ax, 50)
            set(ax,'FontSize',8)
            set(ax,'XLim',[0 25])
            legend(ax, {'True','Predicted'},'FontSize',10)
            legend(ax, 'show');
            ylabel(ax, 'Freq. %', 'FontSize',10);
            title(ax, figTtl, 'FontSize', 10);
            set(fig, 'Visible', 'on');
        end

        function [fig, mdnF1, mnF1, weightedF1]=F1Sizes( ConfusionMat, ...
                figTtl, labels, ignore1st, followed, where, clrs)
            if nargin<7
                clrs=[];
                if nargin<6
                    where=[];
                    if nargin<5
                        followed=[];
                    end
                end
            end
            fig=[];
            N=length(labels);
            if isempty(clrs)
                clrs=zeros(N,3);
                [~,I]=sort(upper(labels));
                for i=1:N
                    idx=I(i);
                    clrs(idx,:)=Gui.HslColor(idx,N);
                end
            end
            if isnumeric(labels)
                if any(labels==0)
                    labels=labels(labels~=0);
                    N=length(labels);
                    if size(ConfusionMat,1)==N+1
                        %argue with argument
                        ignore1st=true;
                    end
                end
            end
            % F1-score
            [f1s, mdnF1, mnF1, sizes, weightedF1]=...
                MatBasics.ConfusionF1(ConfusionMat, ignore1st);
            if isempty(figTtl)
                return;
            end
            totalSize=max(sizes);
            disp(['Median F1-score = ' num2str(mdnF1)])
            fig=Gui.NewFigure(true, false, true);
            X=log10(sizes);
            Y=f1s;
            ax=Gui.Axes(fig);
            hold(ax, 'all');
            hold(ax, 'on');
            Hs=zeros(1,N);
            for i=1:N
                Hs(i)=plot(ax, X(i), Y(i), 'o', 'MarkerFaceColor', ...
                    clrs(i,:), 'MarkerSize', 8);
            end
            %H=scatter(ax, X(l), Y(l), 100, clrs(l,:), 'filled');
            xlim(ax, [0 log10(max(sizes)*1.2)]);
            ylim(ax, [-.05, 1.05]);
            l=cell(1, N);
            txt=['<html><font color="#FF22FF">F1-score</font>' ...
                '/<font color="blue">sizes</font>-subset</html>'];
            [~,I]=sort(f1s, 1, 'descend');
            for k=1:N
                idx=I(k);
                if isnumeric(labels(idx))
                    label=num2str(labels(idx));
                else
                    label=labels{idx};
                end
                l{k}=sprintf(['<html>%s<font color="#FF22FF">%s</font>/' ...
                    '<font color="blue">%s</font>- %s</html>'], ...
                    Html.Symbol(clrs(idx,:), sizes(idx)/totalSize*200, ...
                    false, [10 14]), ...
                    String.encodeRounded(f1s(idx), 3), ...
                    String.encodeK(sizes(idx)), label);
            end
            [~, ~, ~, dlg]=Gui.Ask(struct('msg', txt, ...
                'checkBoxFnc', @(h, e, idx, checkBoxes)flash(idx), 'modal', ...
                false), l, '', '',  0, [], true, N);

            title(ax, {figTtl, sprintf(['F1-score  median/mean ' ...
                '\\fontsize{14}'...
                '\\color{magenta}%s\\color{black}/' ...
                '\\color{magenta}%s\\color{black}'], ...
                String.encodeRounded(mdnF1, 3), ...
                String.encodeRounded(mnF1, 3))});
            xlabel(ax, 'Subset size (log10)');
            ylabel(ax, 'F1-score');
            box(ax, "on");
            grid(ax, "on");
            if ~isempty(followed)
                SuhWindow.Follow(fig, followed, where, true);
                SuhWindow.SetFigVisible(fig);
            else
                set(fig, 'Visible', 'on');
            end
            SuhWindow.Follow(dlg, fig, 'south west+', true);

            function flash(idx)
                idx2=I(idx);
                H=Hs(idx2);
                sz=get(H, 'MarkerSize');
                set(H, 'MarkerSize', sz*2);
                Gui.FlashN(H, 5, .25, false);
                MatBasics.RunLater(@(h,e)resize(H, sz), 2);
            end

            function resize(H, sz)
                set(H, 'MarkerSize', sz);
            end
        end

        function [fig, mdnScore, mnScore, weightedScore, clrs]=ScoreSizes( ...
                scoreTable, scoreColumnIndex, sizeColumnIndex, labels, ...
                figTtl, scoreName, followed, where, clrs, ax, plotMap, args)
            if nargin<12
                args=struct('plotMdnMn', true);
                if nargin<11
                    plotMap=[];
                    if nargin<10
                        ax=[];
                        if nargin<9
                            clrs=[];
                            if nargin<8
                                where=[];
                                if nargin<7
                                    followed=[];
                                    if nargin<6
                                        scoreName='F1-score';
                                    end
                                end
                            end
                        end
                    end
                end
            end
            fig=[]; 
            N=length(labels);
            if isempty(clrs)
                clrs=zeros(N,3);
                [~,I]=sort(upper(labels));
                for i=1:N
                    idx=I(i);
                    clrs(idx,:)=Gui.HslColor(i,N);
                end
            else
                if size(clrs, 1)~=N
                    error('%d colors does not equal %d labels', ...
                        size(clrs, 1), N);
                end
                if iscell(clrs) && startsWith(clrs{1}, '<html>')
                    strs=clrs;
                    clrs=zeros(N,3);
                    for i=1:N
                        if ~isempty(strs{i})
                            clrs(i,:)=Gui.ConvertFromHexColor(...
                                regexp(strs{i}, 'color="#(.*)"', 'tokens'));
                        else
                            clrs(i,:)=[.8 .8 .8];
                        end
                    end
                end
            end
            if isnumeric(labels)
                if any(labels==0)
                    labels=labels(labels~=0);
                    N=length(labels);                    
                end
            end            
            [scores, mdnScore, mnScore, sizes, weightedScore]=...
                MatBasics.ScoresAndSizes(scoreTable, ...
                scoreColumnIndex, sizeColumnIndex);
            if isempty(figTtl)
                return;
            end
            app=BasicMap.Global;
            totalSize=max(sizes);
            disp(['Median ' scoreName '=' num2str(mdnScore)])
            insideOneFig=isempty(ax);
            if insideOneFig
                [fig, tb]=Gui.NewFigure(true, false, true);
                Gui.AddSvgToToolBar(fig, tb);
                set(fig, 'Name', [scoreName ', ' figTtl]);
                ax=Gui.Axes(fig);
            end
            X=log10(sizes);
            Y=scores;
            hold(ax, 'all');
            hold(ax, 'on');
            scoreHs=zeros(1,N);
            if app.highDef
                ms=6;
            else
                ms=8;
            end
            if isempty(plotMap)
                for i=1:N
                    scoreHs(i)=plot(ax, X(i), Y(i), 'o', 'MarkerFaceColor', ...
                        clrs(i,:), 'MarkerSize', ms);
                end
            else
                prefix=[figTtl '.' scoreName '.'];
                for i=1:N
                    scoreHs(i)=plot(ax, X(i), Y(i), 'o', 'MarkerFaceColor', ...
                        clrs(i,:), 'MarkerSize', ms);
                    try
                        plotMap.set([prefix labels{i} '.' num2str(sizes(i))], scoreHs(i));
                    catch
                        plotMap.set([prefix labels{i}], scoreHs(i));
                    end
                end
            end
            %H=scatter(ax, X(l), Y(l), 100, clrs(l,:), 'filled');
            xlim(ax, [0 log10(max(sizes)*1.2)]);
            ylim(ax, [-.05, 1.05]);
            l=cell(1, N);
            txt=['<html><font color="#FF22FF">' scoreName '</font>' ...
                '/<font color="blue">sizes</font>-subset</html>'];
            if insideOneFig
                [~,I]=sort(scores, 1, 'descend');
                for k=1:N
                    idx=I(k);
                    if isnumeric(labels(idx))
                        label=num2str(labels(idx));
                    else
                        label=labels{idx};
                    end
                    html1=Html.EncodeSort('score', scores(idx));
                    html2=Html.EncodeSort('size', sizes(idx));
                    html3=Html.EncodeSort('name', label);
                    l{k}=sprintf(['<html>%s<font color="#FF22FF">%s</font>/' ...
                        '<font color="blue">%s</font>- %s</html>'], ...
                        Html.Symbol(clrs(idx,:), sizes(idx)/totalSize*200, ...
                        false, [10 14]), ...
                        [html1 String.encodeRounded(scores(idx), 3)], ...
                        [html2 String.encodeK(sizes(idx))], ...
                        [html3 label]);
                end
                legendTitle=['Legend: ' strrep(figTtl, ' run on ', ',') ...
                    ',' scoreName];
                [~, ~, ~, dlg]=Gui.Ask(struct('msg', txt, ...
                    'checkBoxFnc', @(h, e, idx, checkBoxes)flash(idx), ...
                    'modal', false), l, '', legendTitle,  ...
                    0, [], true, 8);
            end
            if app.highDef
                fs='\\fontsize{12}';
            else
                fs='\\fontsize{13}';
            end
            if insideOneFig
                if strcmpi(scoreName, 'F1-Score')
                    sn='F1-score ^{AKA F-measure)}';
                else
                    sn=scoreName;
                end
                if ~args.plotMdnMn
                    title(ax, figTtl);
                else
                    ttlMdnMn=sprintf([sn ': median/mean/weighted ' fs...
                        '\\color{magenta}%s\\color{black}/' ...
                        '\\color{magenta}%s\\color{black}/'...
                        '\\color{magenta}%s\\color{black}'], ...
                        String.encodeBank0(mdnScore), ...
                        String.encodeBank0(mnScore),...
                        String.encodeBank0(weightedScore));
                    title(ax, {figTtl, ttlMdnMn});
                end
            else
                rgb=ClassificationTable.RGB;
                if ~args.plotMdnMn
                    title(ax, figTtl);
                else

                    ttlMdnMn={figTtl, sprintf([fs ...
                        rgb '%s\\color{black}/' ...
                        rgb '%s\\color{black}/'...
                        rgb '%s\\color{black}'],...
                        String.encodeBank0(mdnScore), ...
                        String.encodeBank0(mnScore),...
                        String.encodeBank0(weightedScore))};
                    title(ax, ttlMdnMn);
                end
            end
            if insideOneFig
                xlabel(ax, 'Number of cells in population (log10)');
                ylabel(ax, scoreName);
            end
            box(ax, "on");
            grid(ax, "on");
            if insideOneFig
                if ~isempty(followed)
                    SuhWindow.Follow(fig, followed, where, true);
                    SuhWindow.SetFigVisible(fig);
                else
                    set(fig, 'Visible', 'on');
                end
                SuhWindow.Follow(dlg, fig, 'south west+', true);
            end
            
            function flash(idx)
               idx2=I(idx);
               figure(fig);
               H=scoreHs(idx2);
               sz=get(H, 'MarkerSize');
               set(H, 'MarkerSize', sz*2);
               Gui.FlashN(H, 5, .25, false);
               MatBasics.RunLater(@(h,e)resize(H, sz), 2);
            end

            function resize(H, sz)
                set(H, 'MarkerSize', sz);
            end
        end


        function [mdnF1, mnF1, weightedF1, fig1, fig2, fig3]=Deconfuse( ...
                ConfusionMat, figTtl, labels, showBars, follow, where)
            if nargin<6
                where=[];
                if nargin<5
                    follow=[];
                    if nargin<4
                        showBars=false;
                    end
                end
            end
            [fig1,  mdnF1, mnF1, weightedF1]=Plots.F1Sizes( ...
                ConfusionMat, figTtl, labels, false, follow, where);
            if showBars
                fig2=Plots.FreqBarPlot(ConfusionMat, labels, figTtl);
                SuhWindow.Follow(fig2, fig1, 'south west+', true);
            end
            fig3=Plots.FreqScatterPlot(ConfusionMat, labels, figTtl);
            SuhWindow.Follow(fig3, fig1, 'south east+', true);
            figure(fig1);            
        end
    end
end