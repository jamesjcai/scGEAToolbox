classdef MDS< handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(Constant)
        NO_MATCH_CLR=[.99 .99 .99]; 
        CLR_EDGE_NORMAL=[.88 .88 .88];
        CLR_EDGE_SPLIT=[.6 .6 0];        
        FREQ_BY_SIZE=false;
        MATCH_DELIMITER=',\bf ';
        MAX_MATCH_FREQ=.1;
    end
    
    properties(SetAccess=private)
        hslColorCount;
        matchColorCount;
        selectedMarkers;
        heatMap4Browser1;
        heatMap4Browser2;
        edges;
        nameProps=[];
        doubleClicker;
        allHs;
        allDetailsHs;
        legendHs;
        distances;
        eigvals;
        sizes;
        longestName=0;
        measurements;%e.g. lattitude and longtitude
        rawMeasurements;
        faceColorDimension=0; 
        edgeColorDimension=0;
        textFaceColorDimension=0; 
        textEdgeColorDimension=0;        
        symbolDimension=0;
        fig
        ax;
        axDetails;
        tb;
        sizeOffsets;
        mouseEar;
        nDimensions;
        dimensionName;
        dimensionItemName;
        nThings;
        thingName;
        showMeasurements;
        showMeasurementHs;
        allMeasurementNames;
        allMeasurements;
        allRawMeasurements;
        measurementNames;
        measurementIdxs;
        fcsIdxs;
        fcsIdxsStr;
        app=BasicMap.Global;
        yRange=0;
        xRange=0;
        xl;
        yl;
        textPolicy=0;
        xyNormalized;
        xyNormalizedC;
        xySubplotIdxs;
        seePositions=0;
        seeAx=0;
        uniqueSeePositions;
        subplotRows;
        subplotCols;
        axes={};
        btn1D;
        dimFace;
        pc;
        nPc;
        mins;
        maxs;
        window;
        doArrows=[];
        pnlSubPlots;
        btnWrench;
        cb0;
        matchSubsetStrs={};
        matchCnts=[];
        jdMrker;
        seeIdxs=[];
        rowClickedCallback;
        columnClickedCallback;
        cellClickedCallback;
        rowHeaderEnteredCallback;
        btnHeatMap;
    end
    
    properties
        fncWdcu;
        lastNear=0;
        lastTipFnc;
        fncSvg;
        fncToolBar;
        fnc1DPathFinder;
        plotIdxs;
        multiProps;
        positionPrefix='L';
        names;
        matchProps;
        matchIdxs;
        matchFreqs;
        isTop;
        posProps;
        colors=[];
        fncRightClick;
        fncLegendClick;
        useComboBoxes=true;
        fncRefresh=[];
        btnGenie;
        jBtnGenie;
        freqs;
        associatedThis;
        fncDrop;
        globalColorsByName;
        XY;
        ignoreScatter=false;
    end
    
    methods(Static)
        function subsetTxt=BfPrefixToSupSuffix(subsetTxt)
            
            idx=String.IndexOf(subsetTxt, '\bf');
            N=length(subsetTxt);
            if idx>0 
                ending='';
                if N>=idx+3
                    start=subsetTxt(idx+3:end);
                    if idx>1
                        ending=[' <sup>' subsetTxt(1:idx-1) '</sup>'];
                    end
                else
                    start='';
                end
                subsetTxt=[start ending];
            end
        end
        function pc=ColorSet
            N=1024;
            pc=jet(N);%this.app.pseudoColor;
            %pc=winter(N);
            %pc=parula(N);
            %pc=cool(N);
            %pc=spring(N);
            %pc=summer(N);
            %pc=autumn(N);
            %pc=bone(N);
            %pc=copper(N);
        end
        function [r,c]=GetSubPlots(N, cols)
            if N>=cols
                c=cols;
            else
                c=N;
            end
            r=floor(N/cols);
            if mod(N, cols)>0
                r=r+1;
            end
        end
    end
    
    methods
        function setSeePositions(this, seePositions)
            this.seePositions=seePositions;
            usp=unique(seePositions);
            [this.subplotRows, this.subplotCols]=MDS.GetSubPlots(...
                length(usp), 3);
        end
        
        function dimension=findDimension(this, dimension)
            if ischar(dimension)
                if ~isempty(this.measurementNames)
                    dimension=StringArray.IndexOf(this.measurementNames, ...
                        dimension);
                else
                    dimension=0;
                end
            end
        end
        
        function setSymbolDimension(this, dimension)
            this.symbolDimension=this.findDimension(dimension);
        end

        function setTextFaceColorDimension(this, dimension)
            this.textFaceColorDimension=this.findDimension(dimension);
        end

        function setTextEdgeColorDimension(this, dimension)
            this.textEdgeColorDimension=this.findDimension(dimension);
        end

        function setFaceColorDimension(this, dimension)
            this.faceColorDimension=this.findDimension(dimension);
        end

        function setEdgeColorDimension(this, dimension)
            this.edgeColorDimension=this.findDimension(dimension);
        end

        function setNames(this, dimensionItemName, dimensionName, thingName)
            this.dimensionName=dimensionName;
            this.thingName=thingName;
            this.dimensionItemName=dimensionItemName;
        end
        
        function setMouseEar(this, ear)
            this.mouseEar=ear;
        end
        
        function this=MDS(names, measurements, sizes, freqs, ...
                rawMeasurements, globalColorsByName)
            this.multiProps=this.app;
            this.posProps=Map;
            if isequal('java.lang.Object[]', class(names))
                names=StringArray.JavaArrayToCell(names);
            end
            [this.nThings, this.nDimensions]=size(measurements);
            N=this.nThings;
            assert(length(names)==N);
            longest=0;
            for i=1:N
                name=strtrim(names{i});
                if length(name)>longest
                    longest=length(name);
                end
                names{i}=name;
            end
            this.allHs=zeros(1, N);
            this.edges=zeros(N, 4);
            this.allDetailsHs=zeros(1, N);
            this.legendHs=zeros(1, N);
            this.longestName=longest;
            this.names=names;
            this.measurements=measurements;
            if nargin<6
                globalColorsByName=[];
                if nargin<5
                    rawMeasurements=measurements;
                end
            end
            this.globalColorsByName=globalColorsByName;
            this.rawMeasurements=rawMeasurements;
            this.distances=MatBasics.PDist2Self(measurements);     
            [this.XY, this.eigvals] = cmdscale(this.distances);
            this.measurementIdxs=1:length(this.measurements);
            if size(this.XY, 2)==1
                exc= MException('MDS:MDS','Cannot compress distances');
                this.XY
                throw(exc);
            end
            if nargin>3
                [R,C]=size(freqs);
                if R>C
                    this.freqs=freqs';
                else
                    this.freqs=freqs;
                end
            else
                this.freqs=sizes/sum(sizes);
            end
            if nargin>2 && ~isempty(sizes)
                [R,C]=size(sizes);
                if R>C
                    this.sizes=sizes';
                else
                    this.sizes=sizes;
                end
            else
                R=size(this.XY, 1);
                this.sizes=zeros(1,R);
                this.sizes=this.sizes+25;
            end
            this.xyNormalizedC={};
            this.xySubplotIdxs={};
            this.allMeasurements=this.measurements;
            this.allRawMeasurements=this.rawMeasurements;
            this.pc=MDS.ColorSet;
            this.nPc=size(this.pc, 1);
            
        end
        
        function faceColor=getFaceColor(this, idx)
            if ~isempty(this.matchProps) && ~isempty(this.matchIdxs)
                key=[num2str(this.matchIdxs(idx)) '.color'];
                if this.matchProps.has(key)
                    faceColor=str2num(this.matchProps.get(key)); %#ok<ST2NM> 
                    this.matchColorCount=this.matchColorCount+1;
                else
                    mx=max(this.matchIdxs);
                    mi=this.matchIdxs(idx);
                    if mi==0
                        faceColor=MDS.NO_MATCH_CLR;
                    else
                        if ~isempty(this.globalColorsByName)
                            faceColor=this.globalColorsByName.get(...
                                this.names{idx}, mi,mx);
                        else
                            faceColor=Gui.HslColor(mi, mx);
                        end
                        this.hslColorCount=this.hslColorCount+1;
                    end
                end
                return;
            end
            if this.dimFace>0
                idxFace=this.rank(this.dimFace, idx, this.nPc);
                faceColor=this.pc(idxFace, :);
            elseif ~isempty(this.faceColors)
                faceColor=this.faceColors(idx,:);
            else
                if isempty(this.nameProps)
                    this.nameProps=Map;
                end
                name=this.names{idx};
                if this.nameProps.has(name)
                    faceColor=this.nameProps.get(name);
                else
                    n1=size(this.measurements, 1);
                    if ~isempty(this.globalColorsByName)
                        faceColor=this.globalColorsByName.get(...
                            this.names{idx}, idx, n1);
                    else
                        faceColor=Gui.HslColor(idx, n1);
                    end
                    this.nameProps.set(name, faceColor);
                    return;
                end
            end
        end
    end
    
    properties(SetAccess=private)
        faceColors;
        propSeeText;
    end

    methods
        function setSeeTextProperty(this, prop)
            this.propSeeText=prop;
        end
        
        function setColors(this, faceColors)
            if length(this.names)~=length(faceColors) 
                warning('Colors and names must be same length');
            else
                this.faceColors=faceColors;
            end
        end

        function [difOut, stmt]=getMatchFreqDif(this, mi, verbose)
            if mi<=length(this.matchFreqs)
                m1=this.matchFreqs(mi+1,1);
                m2=this.matchFreqs(mi+1,2);
                mx=max([m1 m2]);
                dif_=abs(m1-m2);
                dif=dif_/mx;
                mxf=MDS.MAX_MATCH_FREQ;
                if mx<.05
                    mxf=.37;
                elseif mx<.1
                    mxf=.25;
                elseif mx<.3
                    mxf=.175;
                end
                if dif<mxf
                    difOut=0;
                else
                    difOut=dif;
                end
                    
                if nargout>1
                    if nargin<3
                        if dif>mxf
                            start='<font color="red">';
                            stop='</font>';
                        else
                            start='';
                            stop='';
                        end
                        stmt=sprintf('%s %s relative dif %s between %s & %s', ...
                            start,...
                            String.encodePercent(dif_, mx, 1), ...
                            stop,...
                            [String.encodeRounded(m1/1*100,3,true) '%'],...
                            [String.encodeRounded(m2/1*100,3,true) '%']);
                    else
                        if verbose==1
                            vs=[String.encodeRounded(m2/1*100,1,true) '%'];
                        else
                            vs=[String.encodeRounded(m1/1*100,1,true) '%'];
                        end
                        stmt=sprintf(', ^{vs %s}', vs);
                    end
                end
            else
                stmt='';
                difOut=0;
            end
        end
        function initMultiplePlots(this)
            this.xyNormalizedC={};
            this.xySubplotIdxs={};
            this.axes={};
            this.allDetailsHs=zeros(1, this.nThings);
            this.legendHs=zeros(1, this.nThings);
        end
        
        function setMeasurementNames(this, names, filtered, fcs)
            this.measurementNames=names;
            this.allMeasurementNames=this.measurementNames;
            if nargin>2 && ~isempty(filtered)
                if nargin>3
                    filtered=StringArray.Sort(filtered, fcs.statisticParamNames);
                    this.fcsIdxs=fcs.findFcsIdxs(filtered);
                    this.fcsIdxsStr=MatBasics.toString(this.fcsIdxs,'-');
                end
                this.redoDistances(filtered);
            elseif nargin>3
                this.fcsIdxs=fcs.findFcsIdxs(names);
            end
            this.setComboBoxes(this.fig);
        end
        
        function l=getAxisLabel(this, includeDimensionName)
            if includeDimensionName
                l=this.dimensionName;
            else
                l='';
            end
            if this.nDimensions>2
                l=sprintf('%s %dD\\rightarrow2D', l, this.nDimensions);
            end
        end
        
        function ls=getAxisLabels(this, includeDimensionName)
            axisLabel_ =this.getAxisLabel(includeDimensionName);
            if this.textPolicy==0
                word=' \rm ^{symbol size=freq}';
            else
                word='';
            end
            ls={['\bfMDS1 \rm' axisLabel_  word], ...
                    ['\bfMDS2 \rm' axisLabel_]};
        end

        
        function str=subsetSymbol(this, idx, ~)
            str='';
            if this.allHs(idx)~=0 && ishandle(this.allHs(idx))
                c=get(this.allHs(idx), 'MarkerFaceColor');
                str=['<font  ' Gui.HtmlHexColor(c)...
                    '>&bull;</font>'];
                f=get(this.allHs(idx), 'markerSize');
                f=((f-7)/4)+3;
                if f>8
                    f=8;
                end
                f=String.encodeInteger(ceil(  f) );
                str=['<font size="' f '">' str '</font>'];
            end
        end
        
        function [fldr, matchName]=getLocalImageFldr(this)
            matching=~isempty(this.matchProps);
            if matching
                [fldr, matchName, ~]=fileparts(this.matchProps.propertyFile);
                fldr=fileparts(fldr);
                try
                    fldr=this.app.stripIfUnderUrlFolder(fldr, CytoGate.Folder);
                catch
                end
                fldr=fullfile(fldr, 'matches');
                File.mkDir(fldr);
                fldr=fullfile(fldr, 'images');
                %File.mkDir(fldr);
            else 
                fldr=[];
            end            
        end
        
        
        %used by heat map
        function rowClicked(this, ax, row)
            feval(this.mouseEar, ax, -1, row, nan, nan);
        end
        
        %used by heat map
        function cellClicked(this, ~, row, ~)
            feval(this.fnc1DPathFinder, row, 'south west+');
        end

        %used by heat map
        function columnClicked(this, ax, column)
            fprintf('No SuhHeatMap column listener column=%d/%d, ax=%d\n',...
                column, length(this.measurementIdxs), ~isempty(ax));
        end

        %used by heat map
        function [tip, botPnl]=rowHeaderEntered(this, ax, row)
            [tip, botPnl]=feval(this.mouseEar, ax, true, row, nan, nan);
        end
        
        function [html_, htmlFile]=showMarkerHeatMap(this, doHtmlOnly)
            try
                html_=[];
                htmlFile=[];
                matching=~isempty(this.matchProps);
                if matching
                    showNoMatch=this.matchProps.is('showNoMatch', true);
                else
                    showNoMatch=true;
                end
                if nargin>1 && doHtmlOnly
                    html1=[];
                    html2=[];
                    tp1=[];
                    tp2=[];
                    [html_, htmlFile]=browse(false);
                    return;
                end
                html1=[]; html2=[];
                if ~isempty(this.jdMrker)
                    this.jdMrker.setVisible(true);
                    this.jdMrker.toFront;
                    return;
                end
                if isempty(this.showMeasurements)
                    this.showMeasurements=false(1, length(this.measurementNames));
                end
                sm1=this.app.smallStart;
                sm2=this.app.smallEnd;
                [~,C]=size(this.measurements);
                options=cell(1,C);
                for col=1:C
                    options{col}=['<html><font color="blue"><i>' ...
                        num2str(col) '.</i></font> ' ...
                        this.measurementNames{col}  '</html>'];
                end
                heatMap=Gui.BorderPanel(10,1);
                [top, ~, cb]=SuhHeatMap.DropDowns(@(h,e)setMap(false), ...
                    'mdsHeatMapOrder', @(h,e)setMap(false), 4,{}, ...
                    'East', 'West', ' low to high expression');                
                heatMap.add(top, 'North');
                jd=[];
                tp1=[];
                tp2=[];
                maxN=0;                
                setMap(false);
                nItems=12;
                if maxN<7
                    nItems=9;
                elseif maxN<12
                elseif maxN<20
                    nItems=20;
                else
                    nItems=25;
                end
                this.rowClickedCallback=@(ax, row)...
                    rowClicked(this,ax,row);
                this.columnClickedCallback=@(ax,column)...
                    columnClicked(this,ax,column);
                this.cellClickedCallback=@(ax,row,column)...
                    cellClicked(this, ax, row, column);
                this.rowHeaderEnteredCallback=@(ax,row)...
                    rowHeaderEntered(this, ax, row);
                browseBtn=Gui.NewBtn('Browse', @(h,e)browse(true),...
                    'Click to see in default browser', 'world_16.png');
                sw=Gui.Panel;
                sw.add(browseBtn);
                sw=Gui.AddTipImgCheckBox(this.app, sw, this.lastTipFnc, ...
                    'East', 'West');         
                [opts,cdI]=Fcs.SortXml(options);
                [~,~,jd]=mnuMultiDlg(struct(...
                    'allMsg', ['<html><b>All</b><br>' sm1 ...
                    '(Selections show in<br> MDS labels)' sm2 ...
                    '<hr></html>'],...
                    'noCancel', true, 'where', 'west++',...
                    'icon', 'none', 'modal', false, ...
                    'msg', '', 'checkBoxFnc', ...
                    @(h, e, idxs, checkBoxes)markerSelected(idxs)), ...
                    'Median marker exppngression heat map', opts, ...
                    find(this.showMeasurements)-1, false, true, sw, ...
                    'south west buttons', heatMap, 'West', nItems);
                this.jdMrker=jd;
            catch ex
                BasicMap.Global.reportProblem(ex);
            end
            
            function markerSelected(idxs)
                this.updateMrkrExpr(cdI(idxs));
                setMap(false);
            end
            function [html, htmlFile]=browse(show)
                htmlFile=[];
                if ~show
                    [imgFldr, pngFile]=this.getLocalImageFldr;
                    [imgParentFldr, imgFldr2]=fileparts(imgFldr);
                    htmlFile=fullfile(imgParentFldr, [pngFile '.html']);
                    pngFile=[pngFile '.png'];
                    setMap(true, 3, imgFldr);                    
                else
                    setMap(true);
                    pngFile='lastMds.png';
                    imgFldr=this.app.appFolder;
                    imgFldr2=imgFldr;
                end
                pngPath=fullfile(imgFldr, pngFile);
                saveas(this.fig, pngPath);
                html=['<table><tr><td colspan="2"><br><hr>' ...
                    Html.ImgForBrowser(pngFile, imgFldr2)...
                    '</td></tr></table><table cellpadding="14">'...
                    '<tr><td valign="top">' html1 '</td><td valign="top">' ...
                    html2 '</td></tr></table>'];
                if show
                    Html.Browse(['<html>' html '</html>']);
                end
            end
            
            function setMap(forBrowser, idx, imgFldr)
                if nargin<2
                    idx=cb.getSelectedIndex+1;
                end
                subsetOrder=MarkerHeatMap.ORDER{idx};
                if ~isempty(tp1)
                    if ~forBrowser
                        heatMap.remove(tp1);
                    end
                end
                [I, ax_, are2plots]=this.sortXyAndAx(1, subsetOrder);
                maxN=length(I);
                if nargin>2
                    [html1, tp1]=SuhHeatMap.Html(this, ax_, ...
                        showNoMatch, I, [], true, forBrowser, ...
                        true, 'Marker expression low to high', '0', ...
                        '0', '0', [], [], imgFldr);
                else
                    [html1, tp1]=SuhHeatMap.Html(this, ax_,  ...
                        showNoMatch, I, [], true, forBrowser);
                end
                if are2plots
                    if ~forBrowser
                        heatMap.add(tp1, 'West');
                    end
                    if ~isempty(tp2)
                        if ~forBrowser
                            heatMap.remove(tp2);
                        end
                    end
                    [I, ax_]=this.sortXyAndAx(2, subsetOrder);
                    if nargin>2
                        [html2, tp2]=SuhHeatMap.Html(this, ax_, ...
                            showNoMatch, I, [], true, forBrowser,...
                            true, 'Marker expression low to high', '0', ...
                            '0', '0', [], [], imgFldr);
                    else
                        [html2, tp2]=SuhHeatMap.Html(this, ax_, ...
                            showNoMatch, I, [], true, forBrowser);
                    end
%                    HH=MarkerHeatMap.HtmlOld(this, ax_, showNoMatch, I, [],...
%                        true, forBrowser);                        
                    if ~forBrowser
                        heatMap.add(tp2, 'East');
                    end
                    if length(I)>maxN
                        maxN=length(I);
                    end
                    if forBrowser
                        this.heatMap4Browser2=html2;
                    end
                elseif ~forBrowser
                    this.heatMap4Browser1=html1;
                    heatMap.add(tp1, 'Center');
                end
                if ~forBrowser
                    if ~isempty(jd)
                        SuhHeatMap.Repack(jd, false);
                    end
                else
                    this.heatMap4Browser1=html1;
                end
            end
        end

        function [I, ax_, are2plots]=sortXyAndAx(this, plotIdx, distFromWhere)
            if isempty(this.axes)
                ax_=this.ax;
            else
                ax_=this.axes{plotIdx};
            end
            distFromWhere=char(edu.stanford.facs.swing.Basics.RemoveXml(...
                distFromWhere));
            are2plots=length(this.xyNormalizedC)==2;
            if are2plots
                xy=this.xyNormalizedC{plotIdx};
                rowIdxs=this.plotIdxs==plotIdx;
            else
                xy=this.xyNormalized;
                rowIdxs=1:length(this.names);
            end
            if startsWith(distFromWhere, 'Ascending')
                I=sort_(this.names(rowIdxs));
            elseif startsWith(distFromWhere, 'Descending')
                I=sort_(this.names(rowIdxs));
                I=flip(I);
            elseif startsWith(distFromWhere, 'Highest') ||...
                    startsWith(distFromWhere, 'Lowest')
                I=SuhHeatMap.SortByMeasurement(distFromWhere, this, rowIdxs);
            else
                I=SuhHeatMap.SortDists(ax_, xy, distFromWhere);
            end
            rowIdxs=find(rowIdxs);
            I=rowIdxs(I);
            this.names(I)
            
            function I=sort_(strs)
                out=StringArray.ToLower(strs);
                N=length(out);
                for i=1:N
                    out{i}=strrep(out{i}, '\bf', '');
                end
                [~, I]=sort(out);
            end
        end
        
        function updateMrkrExpr(this, idxs)
            J=jet;
            nJ=size(J,1);
            matching=~isempty(this.matchProps);
            if matching
                showNoMatch=this.matchProps.is('showNoMatch', true);
            else
                showNoMatch=true;
            end
            if nargin>1
                this.selectedMarkers=idxs;
            else
                idxs=this.selectedMarkers;
            end
            if any(idxs)
                this.showMeasurements(:)=false;
                this.showMeasurements(idxs)=true;
            end
            [R,C]=size(this.measurements);
            if isempty(this.showMeasurementHs)
                this.showMeasurementHs=zeros(1, R);
            end
            if ~any(this.showMeasurements)
                for row=1:R
                    if this.showMeasurementHs(row)~=0
                        delete(this.showMeasurementHs(row));
                        this.showMeasurementHs(row)=0;
                    end
                end
                return;
            end
            bgc=[.74 .74 .74];
            for row=1:R
                if ~showNoMatch
                    if this.matchIdxs(row)==0
                        continue;
                    end
                end
                strs={''};
                done=0;
                for col=1:C
                    if ~this.showMeasurements(col)
                        continue;
                    end
                    done=done+1;
                    pk=this.measurements(row,col);
                    cIdx=floor(pk*nJ);
                    if cIdx<1
                        cIdx=1;
                    end
                    c=J(cIdx,:);
                    if mod(done, 7)==0
                        strs{end+1}='';
                    end
                    strs{end}=[strs{end} '\color[rgb]{' num2str(c(1)) ...
                        ',' num2str(c(2)), ',' num2str(c(3)) '}' ...
                        num2str(col) ' '];
                end
                if this.showMeasurementHs(row)==0 ...
                        || ~ishandle(this.showMeasurementHs(row))
                    if isempty(this.axes)
                        ax_=this.ax;
                    else
                        ax_=this.axes{this.seePositions(row)};
                    end
                    [nudgeX, nudgeY]=SubsetLabel.Nudge(ax_, .06, .04);
                    tX=this.XY(row,1)+nudgeX;
                    tY=this.XY(row,2)-nudgeY;
                    mfc=get(this.allHs(row), 'markerFaceColor');
                    this.showMeasurementHs(row)=text(ax_, tX, tY, strs, ...
                        'lineWidth', 3, ...
                        'backgroundColor', bgc, 'edgeColor', mfc);
                    draggable(this.showMeasurementHs(row));
                else
                    set(this.showMeasurementHs(row), 'String', strs);
                end
            end
        end
        
        function [H, ax_]=see(this, seePosition, axisLabels, titleLabel)
            if nargin<2
                seePosition=0;
            end
            this.hslColorCount=0;
            this.matchColorCount=0;
            ax_=this.ax;
            nThings_=this.nThings;
            if length(this.seePositions)>1
                if this.subplotRows>0 && this.subplotCols>0
                    nThings_=sum(this.seePositions==seePosition);
                    p=seePosition;
                    if ~isempty(p)
                        m=this.subplotRows;
                        n=this.subplotCols;
                        if isempty(this.pnlSubPlots)
                            this.pnlSubPlots=uipanel('Parent', this.fig, ...
                                'Units', 'normalized',...
                                'BackgroundColor', 'white',...
                                'ShadowColor', 'white',...
                                'Position', [.03 .04 .955, .945]);
                        end
                        [X,Y]=ind2sub([this.subplotCols this.subplotRows ], p);
                        ax_=subplot(m,n, p, 'Parent', this.pnlSubPlots);
                        Gui.ShrinkNormalized(ax_, 0, .18);
                        margin=.12;
                        firstMargin=.11;
                        margins=firstMargin+(margin*this.subplotRows);
                        sz=(1-margins);
                        height=sz/this.subplotRows;
                        y=firstMargin+(margin*(Y-1))+((Y-1)*height);
                        firstMargin=.06;
                        margin=.055;
                        margins=firstMargin+(margin*this.subplotCols);
                        sz=(1-margins);
                        width=sz/this.subplotCols;
                        x=firstMargin+(margin*(X-1))+( (X-1)*width);
                        if this.subplotCols > this.subplotRows 
                            height=height/this.subplotRows;
                        elseif this.subplotRows>this.subplotCols
                            width=width/this.subplotCols;
                        end
                        %P=get(ax_, 'position');
                        set(ax_, 'Position', [x y width height], 'Parent', this.pnlSubPlots, 'units', 'normalized');
                        %ax_=subplot(m,n, p, 'Parent', this.fig, 'Position', [P(1) P(2) sz sz]);
                        %ax_=subplot('Position', [P(1) P(2) sz sz], 'Parent', this.fig);
                        get(ax_, 'Position')
                        this.axes{end+1}=ax_;
                    end
                end
            end
            if nargin<4
                hasDimName=~isempty(this.dimensionName);
                if ismac
                    fs1='\fontsize{14}';
                    fs2='\fontsize{11}';
                else
                    fs1='\fontsize{14}';
                    fs2='\fontsize{9}';
                end
                if ~isempty(this.thingName)
                    if hasDimName
                        titleLabel={sprintf('%s%d %s', fs1, ...
                            nThings_, this.thingName), ...
                            sprintf('%sin %dD space (%s)', fs2, ...
                            this.nDimensions, this.dimensionName)};
                    else
                        titleLabel={sprintf('%s%d %s', fs1, ...
                            nThings_, this.thingName), ...
                            sprintf('%sin %dD space', fs2, ...
                            this.nDimensions)};
                    end
                elseif hasDimName
                    titleLabel=sprintf('%s%dD space (%s)', ...
                            fs1, this.nDimensions, this.dimensionName);
                else
                    titleLabel=[];
                end
                if nargin<3
                    axisLabels=this.getAxisLabels(true);
                end
            end
            cla(ax_, 'reset');
            
            hold(ax_, 'all');            
            P=Gui.GetPixels(ax_);
            pixelWidth=P(3);
            MIN=60; % 5 pixels ( sqrt(100)==10 )
            MAX=1500; % 40 pixels
            mn=min(this.sizes);
            sizeRange=range(this.sizes);
            if sizeRange==0
                sizeRange=1;
            end
            fs=12;
            if this.longestName>12
                fs=10;
            end
            if MDS.FREQ_BY_SIZE
                frqs=(this.sizes-mn)/sizeRange;
            else
                frqs=this.freqs/max(this.freqs);
            end
            symbolSize=floor(MIN+(frqs*MAX));
            symbolSize=sqrt(symbolSize);
            maxSymbolSize=max(symbolSize);
            legendSymbolSize=floor((MIN/2.5)+(frqs*(MAX/2.5)));
            legendSymbolSize=sqrt(legendSymbolSize);
            Hs=[];
            legendHs_=[];
            dimSym=getDim(this.symbolDimension);
            this.dimFace=getDim(this.faceColorDimension);
            dimEdge=getDim(this.edgeColorDimension);
            dimTextFace=getDim(this.textFaceColorDimension);
            dimTextEdge=getDim(this.textEdgeColorDimension);
            this.maxs=max(this.measurements);
            this.mins=min(this.measurements);
            symbols={'o', 's', 'd', 'h', '^', 'v', '>', '<', 'p', '*', '+'};
            nSyms=length(symbols);
            N=this.nThings;
            totalSize=0;
            matching=~isempty(this.matchProps);
            if matching
                showNoMatch=this.matchProps.is('showNoMatch', true);
            else
                showNoMatch=true;
            end
            for i=1:N
                totalSize=totalSize+this.sizes(i);
                if length(this.seePositions) >= i && ...
                        this.seePositions(i) ~= seePosition
                    continue;
                end
                if ~showNoMatch
                    if this.matchIdxs(i)==0
                        continue;
                    end
                end
                if dimSym>0
                    idxSym=this.rank(dimSym, i, nSyms);
                else
                    idxSym=1;
                end
                faceColor=this.getFaceColor(i);
                if dimEdge>0
                    idxEdge=this.rank(dimEdge, i, this.nPc);
                    edgeClr=this.pc(idxEdge, :);
                else
                    if matching && ~isempty(this.matchIdxs)
                        edgeClr=[];
                        mi=this.matchIdxs(i);
                        key=[num2str(mi) '.madUnits'];
                        if this.matchProps.has(key)
                            madUnits=str2num(this.matchProps.get(key)); %#ok<ST2NM> 
                            madUnits=madUnits(~isnan(madUnits));
                            if any(madUnits>3)
                                idxSym=4;
                            elseif any(madUnits>2)
                                idxSym=9;
                            end
                            nMatches=this.matchProps.getNumeric([num2str(mi) '.count'],0);
                            if nMatches>2
                                edgeClr=MDS.CLR_EDGE_SPLIT;
                            end
                        else
                            idxSym=2;
                            edgeClr=[1 0 0];
                        end
                        if false %avoid square idxSym==1
                            dif=this.getMatchFreqDif(mi);
                            if dif>0
                                idxSym=2;
                            end
                        end
                    else
                        edgeClr=[.05 .05 .05];
                    end
                end
                if isempty(edgeClr)
                    edgeClr=MDS.CLR_EDGE_NORMAL;
                    lineWidth=1;
                else
                    if symbolSize(i)<maxSymbolSize/4
                        lineWidth=1;
                    elseif symbolSize(i)<maxSymbolSize/2
                        lineWidth=2;
                    else
                        lineWidth=3;
                    end

                end
                if this.textPolicy==1 || this.textPolicy==4
                    legendHs_(end+1)=plot(ax_, this.XY(i, 1), ...
                        this.XY(i, 2), symbols{idxSym}, ...
                        'markerSize', legendSymbolSize(i),...
                        'markerFaceColor', faceColor, ...
                        'markerEdgeColor', edgeClr, 'linewidth', 1);
                    this.legendHs(i)=legendHs_(end);
                end
                H=plot(ax_, this.XY(i, 1), this.XY(i, 2),...
                    symbols{idxSym}, 'markerSize', symbolSize(i),...
                    'markerFaceColor', faceColor, ...
                    'markerEdgeColor', edgeClr, ...
                    'linewidth', lineWidth);
                this.edges(i,:)=[edgeClr lineWidth];
                Hs(end+1)=H;
                this.allHs(i)=H;
            end
            this.seeIdxs=find(this.allHs~=0 & ishandle(this.allHs));
            MX=max(this.XY(:,1:2));
            MN=min(this.XY(:, 1:2));
            if MN(1)<MN(2)
                mnXy=MN(1);
            else
                mnXy=MN(2);
            end
            if MX(1)>MX(2)
                mxXy=MX(1);
            else
                mxXy=MX(2);
            end
            rXy=mxXy-mnXy;
            XL=[mnXy-(.1*rXy) mxXy+(.15*rXy)];
            YL=[mnXy-(.1*rXy) mxXy+(.15*rXy)];
            set(ax_, 'xlim', XL);
            set(ax_, 'ylim', YL);
            xr=.2*(XL(2)-XL(1));
            yr=.2*(YL(2)-YL(1));
            this.initMatchSubsetStrs;
            if this.textPolicy==2%|| dimTextFace>0 || dimTextEdge>0
                nudgeX=SubsetLabel.Nudge(ax_, .07, .02);
                arrows=cell(1, N);
                offX=symbolSize/2/pixelWidth*(XL(2)-XL(1));
                for i=1:N
                    if length(this.seePositions) >= i && ...
                            this.seePositions(i) ~= seePosition
                        continue;
                    end
                    if ~showNoMatch
                        if this.matchIdxs(i)==0
                            continue;
                        end
                    end

                    if dimTextFace>0
                        idxFace=this.rank(dimTextFace, i, this.nPc);
                        backC=this.pc(idxFace, :);
                    else
                        backC='none';
                    end
                    if dimTextEdge>0
                        idxEdge=this.rank(dimTextEdge, i, this.nPc);
                    else
                        idxEdge= 0;
                    end
                    if idxEdge<1
                        edgeC='none';
                    else
                        lw=3;
                        edgeC=this.pc(idxEdge, :);
                        if matching && ~isempty(this.matchIdxs)
                            key=[num2str(this.matchIdxs(i)) '.color'];
                            if this.matchProps.has(key)
                                edgeC=str2num(this.matchProps.get(key)); %#ok<ST2NM> 
                            end
                        end
                    end
                    if mean(backC)< .5
                        color='white';
                    else
                        color='black';
                    end
                    name=this.names{i};
                    del=strfind(name, MDS.MATCH_DELIMITER);
                    if ~isempty(del)
                        str={name(1:del(1)), [name(del(1)+1:end) ',']};
                    else
                        str={name, ''};
                    end
                    if MDS.FREQ_BY_SIZE
                        str{2}=[str{2} '\bf' String.encodePercent(...
                            this.sizes(i), totalSize, 1)];
                    else
                        str{2}=[str{2} '\bf' String.encodePercent(...
                            this.freqs(i), 1, 1)];
                    end
                    key=this.getPositionKey(i);
                    if this.posProps.has(key)
                        pos=str2num(this.posProps.get(key)); %#ok<ST2NM> 
                        tX=pos(1);
                        tY=pos(2);
                    else
                        tX=this.XY(i,1)+(offX(i)+nudgeX);
                        tY=this.XY(i,2);
                    end
                    if strcmp('none', edgeC)
                        T=text(ax_, tX, tY, str, 'FontSize', fs, ...
                            'backgroundColor', backC, 'color', color, ...
                            'HorizontalAlignment', 'center');
                    else
                        T=text(ax_, tX, tY, str, 'FontSize', fs, ...
                            'edgeColor', edgeC, 'backgroundColor', backC, ...
                            'color', color, 'linewidth', 3, ...
                            'HorizontalAlignment', 'center');
                    end
                    drawArrow(T, i, ax_);
                    draggable(T, @(h)motionFcn(h, i, ax_), 'endfcn', ...
                        @(H)endFcn(H, ax_, i),...
                        [XL(1)-xr XL(2)+xr YL(1)-yr YL(2)+yr] );
                end
            elseif this.textPolicy==1 || this.textPolicy==4
                nms={};
                allIdxs=[];
                orderByF=this.textPolicy==4;
                fMeasures=zeros(1, length(this.sizes));
                for i=1:N
                    if length(this.seePositions) >= i && ...
                            this.seePositions(i) ~= seePosition
                            continue;
                    end
                    if ~showNoMatch
                        if this.matchIdxs(i)==0
                            continue;
                        end
                    end
                    if MDS.FREQ_BY_SIZE                        
                        nms{end+1}=[this.names{i} '  \bf' ...
                            String.encodePercent(this.sizes(i), totalSize, 1)];
                    else
                        freq=String.encodePercent(this.freqs(i), 1, 1);
                        dlm='  \bf';
                        if ~isempty(this.matchIdxs) 
                            mi=this.matchIdxs(i);
                            hasF=false;
                            key=[num2str(mi) '.fMeasure'];
                            if this.matchProps.has(key)
                                f=this.matchProps.getNumeric(key, -1);
                                fMeasures(i)=f;
                                if f>=0
                                    hasF=true;
                                    freq=['\color[rgb]{.29 .29 .63}F=' ...
                                        String.encodePercent(f, 1, 1) ', ' freq];
                                end
                            end
                            if ~hasF
                                [diff, stmt]=this.getMatchFreqDif(mi, ...
                                    this.plotIdxs(i));
                                if diff>0
                                    dlm='  \bf\color{magenta}\it';
                                    freq=[freq stmt ];
                                end
                            end
                        end
                        nms{end+1}=[this.names{i} dlm freq];
                    end
                    sMatch=this.getMatchSubsetStr(i, true);
                    nms{end}=[nms{end} ', \color{blue}' sMatch];
                    allIdxs(end+1)=i;
                end
                if ~isempty(this.matchIdxs)
                    mnMatchIdx=min(this.matchIdxs);
                    mxMatchIdx=max(this.matchIdxs);
                    lgndHs=[];
                    plotHs=[];
                    nms_={};
                    allIdxs=[];
                    matches=[];
                    for j=mxMatchIdx:-1:mnMatchIdx
                        k=1;
                        for i=1:N
                            if length(this.seePositions) >= i && ...
                                    this.seePositions(i) ~= seePosition
                                continue;
                            end
                            if ~showNoMatch
                                if this.matchIdxs(i)==0
                                    continue;
                                end
                            end
                            if this.matchIdxs(i)==j
                                lgndHs(end+1)=legendHs_(k);
                                plotHs(end+1)=Hs(k);
                                allIdxs(end+1)=i;
                                nms_{end+1}=nms{k};
                                matches(end+1)=j;
                            end
                            k=k+1;
                        end
                    end
                    %                   [HL, ics, pls txts]=legend(ax_, Hs_, nms_, 'Location', 'northeast','AutoUpdate', 'on');
                    if ~orderByF || all(fMeasures<=0)
                        [~,II]=sort(this.sizes(allIdxs), 'descend');
                    else
                        [~,II]=sort(fMeasures(allIdxs), 'descend');
                    end
                    lgndHs=lgndHs(II);
                    HL=legend(ax_, lgndHs, nms_(II),...
                        'Location', 'northeast','AutoUpdate', 'off');
                    HL.ItemHitFcn=@(h,e)hit(ax_, e, lgndHs, allIdxs, ...
                        II, matches);
                else
                    [~,II]=sort(this.sizes(allIdxs), 'descend');
                    legendHs_=legendHs_(II);
                    HL=legend(ax_, legendHs_, nms(II), ...
                        'Location', 'northeast','AutoUpdate', 'off');
                    HL.ItemHitFcn=@(h,e)hit(ax_, e, legendHs_, allIdxs, ...
                        II);
                end
                P=get(HL, 'Position');
                set(HL, 'Position', [P(1) P(2) P(3) P(4)*1.25]);
                P=get(HL, 'pos');
                alter=1.5;
                height=alter*P(4);
                if height>.97
                    height=.97;
                end
                dif=height-P(4);
                if dif>0
                    Y=P(2)-dif;
                    if Y<.02
                        Y=.02;
                    end
                    set(HL, 'position', [P(1)*1.1 Y P(3)*1.1 height]);
                end
            end
            if ~isempty(axisLabels)
                xlabel(ax_, axisLabels{1});
                if isempty(this.seePositions)||seePosition==1||seePosition==0
                    ylabel(ax_, axisLabels{2});
                end
            end
            if ~isempty(this.seePositions)&& seePosition~=1 && seePosition~=0
                set(ax_, 'Yticklabel', []);
            end
            %set(ax_, 'Xticklabel', []);
            if ~isempty(titleLabel)
                TT=title(ax_, titleLabel);
                %plotedit(TT, 'on');
            end
            this.yRange=YL(2)-YL(1);
            this.xRange=XL(2)-XL(1);
            this.xl=XL;
            this.yl=YL;
            %this.xyNormalized=MDS.Normalize(ax_, this.XY(:,[1 2]));
            if length(this.seePositions) >= i 
                this.xyNormalizedC{end+1}=MDS.Normalize(ax_, this.XY(this.seePositions==seePosition,[1 2]));
                this.xySubplotIdxs{end+1}=find(this.seePositions==seePosition);
            else
                this.xyNormalized=MDS.Normalize(ax_, this.XY(:,[1 2]));
            end
            
            function hit(ax_, event, lgndH, allIdxs, sortII, matches)
                idx=find(lgndH==event.Peer, 1);
                idx=sortII(idx);
                hasEditNameFnc=~isempty(this.fncLegendClick);
                if nargin>5 && matches(idx)>0
                    if hasEditNameFnc
                        if ispc
                            MatBasics.DoLater(@(h,e)flashThem, .9);
                        else
                            MatBasics.DoLater(@(h,e)flashThem, .4);
                        end
                        edit
                        flashThem;
                    else
                        flashThem
                        MatBasics.DoLater(@(h,e)showTip, .3);
                    end
                else
                    if hasEditNameFnc
                        MatBasics.DoLater(@(h,e)flashIt, .4);
                        edit
                    else
                        flashIt;
                    end
                    
                end
                
                function flashIt
                    this.flashIt(allIdxs(idx));
                end
                
                function flashThem
                    this.flashMatches(matches(idx));
                end
                function edit
                    feval(this.fncLegendClick, event.Peer, allIdxs(idx));
                    MatBasics.DoLater(@(h,e)showTip, .3);
                end
                
                function showTip 
                    feval(this.mouseEar, ax_, true, allIdxs(idx));
                end
            end
            
            function drawArrow(T, idx, ax_)
                if isempty(this.doArrows)
                    this.doArrows=this.multiProps.is('ndSeeArrow', false);
                end
                if ~this.doArrows
                    return;
                end
                if ~isempty(arrows{idx})
                    delete(arrows{idx});
                end
                PP=Gui.GetPixels(ax_);
                offX_=symbolSize(idx)/2/PP(3)*(XL(2)-XL(1));
                offY_=symbolSize(idx)/2/PP(4)*(YL(2)-YL(1));
                [xu, yu]=circle(this.XY(idx, 1), ...
                    this.XY(idx, 2), offX_, offY_);
                [cluXy, lblXy]=SubsetLabel.ArrowPts(T, [xu;yu]',...
                    0, 0);
                set(0, 'currentFigure', this.fig);
                set(this.fig, 'currentAxes', ax_);
                arrows{idx}=arrowNewGraphics(lblXy, cluXy, 10, ...
                    'EdgeColor', 'black');
            end
            
            function motionFcn(H, idx, ax_)
                this.app.closeToolTip;
                drawArrow(H, idx, ax_);
            end
            
            function endFcn(H, ax, idx)
                disp('ended');
                this.setMouse;
                if ~isempty(this.mouseEar) && ~DoubleClicker.isMultiSelect
                    feval(this.mouseEar, ax, true, idx, 240, -50);
                end
                key_=this.getPositionKey(idx);
                pos_=get(H, 'position');
                this.posProps.set(key_, num2str(pos_));
                this.posProps.save;
            end
            
            function dim=getDim(dim)
                if dim<1||dim>this.nDimensions
                    dim=0;
                end
            end
            
            function [xu, yu]=circle(x, y, rX, rY)
                th=0:pi/50:2*pi;
                xu=rX*cos(th)+x;
                yu=rY*sin(th)+y;                
            end
            
        end
        function str=getMatchSubsetStr(this, idx, nullIf1x1Or0)
           if idx>length(this.matchSubsetStrs)
               str='';
           else
               str=this.matchSubsetStrs{idx};
               if nargin>2&&nullIf1x1Or0
                   if strcmp('1x1', str)
                       str='';
                   elseif String.StartsWith(str, '0x') || ...
                           String.EndsWith(str, 'x0') 
                       str='';
                   end
               end
           end
        end
        
        
        function setMouse(this)
            set(this.fig,'WindowButtonMotionFcn', @(h,e)motion(this));
            %set(this.fig, 'WindowButtonUpFcn', @(obj,evt)wbu(this, obj));
           
            this.doubleClicker=DoubleClicker2(this.fig, []);
            this.doubleClicker.setFnc([], @(h,cp,handles)wbu(this,h), ...
                [], @(h, cp, handles)wdcu(this, cp), []);
            
        end
        
        function wdcu(this, cp)
            if ~isempty(this.fncWdcu)
                feval(this.fncWdcu, this, cp);
            end
        end
        
        function [plot1Idxs,  plot2Idxs]=getPlotIdxs(this, match)
            plot1Idxs=[];
            plot2Idxs=[];
            allMi=this.matchIdxs;
            for j=1:this.nThings
                if allMi(j)==match
                    plotIdx=this.plotIdxs(j);
                    if plotIdx==1
                        plot1Idxs(end+1)=j;
                    else
                        plot2Idxs(end+1)=j;
                    end
                end
            end
        end
        
        function initMatchSubsetStrs(this)
            if ~isempty(this.matchSubsetStrs)
                return;
            end
            u=unique(this.plotIdxs);
            nU=length(u);
            if nU~=2
                return;
            end
            this.matchSubsetStrs=cell(1, this.nThings);
            this.matchCnts=zeros(this.nThings, 2);
            allMi=this.matchIdxs;
            matches=unique(allMi);
            N=length(matches);
            for i=1:N
                match=matches(i);
                plotCnts=zeros(nU,1);
                idxs=[];
                for j=1:this.nThings
                    if allMi(j)==match
                        plotIdx=this.plotIdxs(j);
                        plotCnts(plotIdx)=plotCnts(plotIdx)+1;
                        idxs(end+1)=j;
                    end
                end
                cntS_1=num2str(plotCnts(1));
                cntS_2=num2str(plotCnts(2));
                str=[ cntS_1 'x' cntS_2 ];
                for j=1:length(idxs)
                    idx=idxs(j);
                    this.matchSubsetStrs{idx}=str;
                    this.matchCnts(idx,:)=plotCnts;
                end
            end
        end
        
        function key=getPositionKey(this, i)
            key=[this.positionPrefix '.' num2str(i) '.position'];
        end
        
        function ranking=rank(this, dim, row, maxRank)
            fprintf('%s dim=%d expr=%s ', this.names{row}, dim, num2str(this.measurements(row, dim)));
            ranking=MatBasics.Rank(this.measurements(row, dim), ...
                this.mins(dim), this.maxs(dim), maxRank);
            fprintf(' rank=%d\n', ranking);
        end
        
        function setLine(this, i)
            num=this.edges(i,:);
            edgeColor=num(1:3);
            lineWidth=num(4);
            set(this.allHs(i), 'lineWidth', lineWidth, ...
                'markerEdgeColor', edgeColor);
        end
        
        function flashMatches(this, match)
            this.app.closeToolTip;
            if match==0
                return;
            end
            try
                plotHs_=[];
                isDetail=[];
                N=length(this.matchIdxs);
                for i=1:N
                    if ~ishandle(this.allHs(i)) || this.allHs(i)==0
                        continue;
                    end
                    if this.matchIdxs(i)==match
                        if isequal('on', get(this.allHs(i), 'visible'))
                            plotHs_(end+1)=this.allHs(i);
                            uistack(plotHs_(end), 'top');
                            isDetail(end+1)=false;
                        elseif ishandle(this.allDetailsHs(i)) ...
                                && this.allDetailsHs(i)~=0
                            plotHs_(end+1)=this.allDetailsHs(i);
                            isDetail(end+1)=true;
                        end
                    end
                    this.setLine(i);
                end
                N2=length(plotHs_);
                for id = 1:3 % Repeat 3 times
                    for i=1:N2
                        set(plotHs_(i), 'visible', 'off');
                        if ~isDetail(i)
                            set(plotHs_(i), 'lineWidth', 4, ...
                                'markerEdgeColor', [0 0 0]);
                        end
                    end
                    pause(0.15);
                    for i=1:N2
                        set(plotHs_(i), 'visible', 'on');
                    end
                    pause(0.15);
                end
            catch ex
                ex.getReport
            end
        end
        
        function flashIt(this, idx)
            this.app.closeToolTip;
            if ~ishandle(this.allHs(idx)) || this.allHs(idx)==0
                return;
            end
            N=this.nThings;
            for i=1:N
                this.setLine(i)
            end
            if isequal('on', get(this.allHs(idx), 'visible'))
                plotHs_=this.allHs(idx);
                uistack(plotHs_(end), 'top');
                set(plotHs_, 'lineWidth', 4, 'markerEdgeColor', [0 0 0])
            else
                plotHs_=this.allDetailsHs(idx);
            end
            for id = 1:3 % Repeat 3 times
                set(plotHs_, 'visible', 'off');
                pause(0.15);
                set(plotHs_, 'visible', 'on');
                pause(0.15);
            end
        end
        
        function keyDown(this, src, evd)
            pressedMeta=StringArray.Contains(evd.Modifier, 'command') && isequal('0', evd.Key);
            pressedCtrl=StringArray.Contains(evd.Modifier, 'control') && isequal('control', evd.Key);
            if pressedCtrl||pressedMeta
                this.app.closeToolTip;
            end
        end
        
        function enable(this, yes)
            Gui.setEnabled(this.fig, yes);
            if ~isempty(this.tb)
                this.tb.setEnabled(yes)
            end
        end

        function str=getMatchesString(this, html1tex2plain3)
            if nargin<2
                html1tex2plain3=2;
            end            
            [unmatched, totals]=this.getMatches;
            str=QfHiDM.GetMatchesString(unmatched, totals, html1tex2plain3);
        end
        
        function [unmatched, total]=getMatches(this)
            unmatched=[];
            total=[];
            
            calc(2)
            calc(1);
            
            function calc(plot)
                if any(this.plotIdxs==plot)
                    total(end+1)=sum(this.plotIdxs==plot);
                    unmatched(end+1)=sum(this.plotIdxs==plot ...
                        & this.matchIdxs==0);
                end
            end
        end
    end
    
    methods(Static)
        function [xy, xRange_, yRange_]=Normalize(ax_, xyIn)
            xl_=xlim(ax_);
            yl_=ylim(ax_);            
            yRange_=yl_(2)-yl_(1);
            xRange_=xl_(2)-xl_(1);
            xy=(xyIn-[xl_(1) yl_(1)])./[xRange_ yRange_];
        end
    end
    
    methods
        function xy=normalize(this, ax, xy)
            xy=(xy-[this.xl(1) this.yl(1)])./[this.xRange this.yRange];
        end

        function [x, y, ax_, axesIdx]=getCp(this)
            if isempty(this.axes) % just ONE ax not 2 axes
                cp = get(this.ax,'CurrentPoint');
                ax_=this.ax;
                x=cp(1, 1);
                y=cp(1, 2);     
                axesIdx=0;
            else
                N=length(this.axes);
                for axesIdx=1:N
                    ax_=this.axes{axesIdx};
                    xl_=xlim(ax_);
                    yl_=ylim(ax_);
                    cp=get(ax_,'CurrentPoint');
                    x=cp(1, 1);
                    y=cp(1, 2);
                    if x>=xl_(1) && x<=xl_(2) && y>=yl_(1) && y<=yl_(2)
                        break;
                    end
                end
            end
            if isempty(cp)
                ax_=[];
            end
        end
        
        function motion(this, h, e)
            [x, y, ax_, idx]=this.getCp;
            if isempty(ax_)
                return;
            end
            this.findNearest(ax_, idx, x, y, true);
        end
        
        function wbu(this, h)
            if ~isempty(this.fncDrop)
                objType=class(get(this.fig, 'UserData'));
                if Gui.IsDroppedOn(this, h, objType, this.fncDrop)
                    return;
                end
            end
            [x, y, ax_, plotIdx]=this.getCp;
            if isempty(ax_)
                return;
            end
            this.app.closeToolTip;
            a=get(gcbf, 'SelectionType');
            if strcmp(a, 'alt')
                if ~isempty(this.fncRightClick)
                    [x, y, ax_, plotIdx]=this.getCp;
                    thingIdx=this.findNearest(ax_, plotIdx, x, y, true);
                    feval(this.fncRightClick, this, thingIdx, x, y, ...
                        ax_, plotIdx, h);
                    return;
                end
            end
            this.findNearest(ax_, plotIdx, x, y, false);
        end
        
        function yes=hasDetails(this, idx)
            yes=this.allDetailsHs(idx)~=0;
        end

        function yes=hasLegend(this, idx)
            yes=this.legendHs(idx)~=0;
        end

        function yes=areDetailsShowing(this, idx)
            yes=this.allDetailsHs(idx)~=0 && isequal('on', ...
                get(this.allDetailsHs(idx), 'visible'));
        end
        
        function toggleDetails(this, idx)
            if ~ishandle(this.allHs(idx)) || this.allHs(idx)==0
                return;
            end
            if this.hasDetails(idx)
                if this.areDetailsShowing(idx)
                    set(this.allDetailsHs(idx), 'visible', 'off');
                    if this.hasLegend(idx)
                        set(this.legendHs(idx), 'visible', 'on');
                    end
                    set(this.allHs(idx), 'visible', 'on');
                    uistack(this.allHs(idx), 'top');
                else
                    set(this.allDetailsHs(idx), 'visible', 'on');
                    set(this.allHs(idx), 'visible', 'off');
                    if this.hasLegend(idx)
                        set(this.legendHs(idx), 'visible', 'off');
                    end
                    uistack(this.allDetailsHs(idx), 'top');
                end
            end
        end
        
        function ax_=getAx(this, idx)
            if this.seePositions==0
                ax_=this.ax;
            else
                ax_=this.axes{this.seePositions(idx)};
            end
        end
        
        function setDetails(this, idx, X, Y, visible)
            if ~ishandle(this.allHs(idx)) || this.allHs(idx)==0
                return;
            end
            if nargin<5
                visible='on';
            end
            ax_=this.getAx(idx);
            mfc=get(this.allHs(idx), 'markerFaceColor');
            if min(mfc)>.8
                mfc=mfc-.2;
            end
            totalColor=sum(mfc);
            if totalColor>1.6
                mfc=mfc*.75; %darken for single cell display                
            end
            if this.allDetailsHs(idx)~=0 && ishandle(this.allDetailsHs(idx))
                delete(this.allDetailsHs(idx));
            end
            try
                if ~isempty(X) && ~isempty(Y)
                    this.allDetailsHs(idx)=plot(ax_, X, Y, '.', ...
                        'markerSize', 1, 'markerEdgeColor', mfc, ...
                        'markerFaceColor', mfc, 'lineStyle', 'none',...
                        'visible', visible);
                end
            catch ex
                ex.getReport
            end
            if isequal(visible, 'on')
                set(this.allHs(idx), 'visible', 'off');
                if this.hasLegend(idx)
                    set(this.legendHs(idx), 'visible', 'off');
                end
            end
            %uistack(this.allDetailsHs(idx), 'bottom');
        end
        
        function setColor(this, idx, clr)
            if ~ishandle(this.allHs(idx)) || this.allHs(idx)==0
                return;
            end
            if this.hasDetails(idx)
                set(this.allDetailsHs(idx), 'markerEdgeColor', clr, ...
                    'markerFaceColor', clr);
            end
            set(this.allHs(idx), 'markerFaceColor', clr);
            if this.hasLegend(idx)
                set(this.legendHs(idx), 'markerFaceColor', clr);
            end
            
        end
    end
    
    methods(Static)
        function ok=DEBUG1
            ok=false;
        end

        function [mds, fig]=New(names, measurementNames, ...
                clrs, mdns, sizes, ttl, locate)
            if nargin<3
                locate=[];
            end            
            mds=MDS(names, mdns, sizes);
            if isempty(locate)
                locate={gcf, 'east', false};
            end
            mds.setColors(clrs);
            [fig, ax]=mds.newFig(ttl, false, locate, ...
                'QfTable.MDS.fig', [], true);
            orlovaScale=.46;
            busy=Gui.ShowBusy(fig, Gui.YellowH3(...
                ['<br>Opening Darya''s tools ...<br>'...
                '(MDS, QFMatch, QF-tree etc.)']),...
                'orlova.png', orlovaScale, false);
            [~, jBtnGenie]=Gui.ImageLabel([], 'smallGenie.png', ...
                '', [], fig, 4, 3);
            Gui.SetTransparent(jBtnGenie);
            mds.setSeeTextProperty('QfTable.MDS.text');
            mds.setMouseEar(@notifyMouse);
            mds.setNames('expression', 'marker', 'Subsets');
            mds.setMeasurementNames(measurementNames);
            mds.see;
            Gui.HideBusy(fig, busy, true);
            app=BasicMap.Global;
            
            function notifyMouse(~, isMotion, idx,x, y)
                if ~isMotion
                    msg(['Clicked on ' ...
                        String.RemoveTex(names{idx})]);
                elseif idx>0
                    if app.highDef
                        x=-40;
                        y=-21;
                        y=y-(33*app.toolBarFactor);
                    else
                        if isempty(ax)
                            x=-51;
                        else
                            x=-31;
                        end
                        y=1;
                    end
                    tip=String.RemoveTex( names{idx} );
                    tipGui.ShowToolTipHere(tip, jBtnGenie,...
                        ax, fig, app, 6, ...
                        [], false, x, y);
                end
            end
        end
    end
    
    methods
        function idx=findNearest(this, ax_, plotIdx, x, y, isMotion, runFnc)
            try
                if nargin<7
                    runFnc=true;
                end
                [xy, xRange_, yRange_]=MDS.Normalize(ax_, [x y]);
                if plotIdx==0
                    near2=pdist2(this.xyNormalized, xy);
                    [minDist,idx]=min(near2);
                else
                    near2=pdist2(this.xyNormalizedC{plotIdx}, xy);
                    [minDist,idx]=min(near2);
                    idxs=this.xySubplotIdxs{plotIdx};
                    idx=idxs(idx);
                    if xy(1)>.96 && xy(1)<.98 && xy(2)>.48 && xy(2)<.50
                        MDS.Normalize(ax_, [x y]);
                    end
                    if MDS.DEBUG1
                        fprintf('xy=%s/%s range=%s/%s\n', ...
                            String.encodeRounded(xy(1), 2, true), ...
                            String.encodeRounded(xy(2), 2, true), ...
                            String.encodeRounded(xRange_, 2, true), ...
                            String.encodeRounded(yRange_, 2, true));
                    end
                end
                if minDist>.1
                    idx=0;
                end
                rx=round(x);
                ry=round(y);
                if idx>0 && idx<=length(this.names)
                    X=this.XY(idx, 1);
                    Y=this.XY(idx, 2);
                    if ~isempty(this.mouseEar) && runFnc
                        if ~isMotion || this.lastNear ~= idx
                            %if this.lastNear ~= idx && isMotion
                            %    fprintf('idx=%d, lastNear=%d\n', idx,this.lastNear);
                            %end
                            this.lastNear=idx;
                            feval(this.mouseEar, ax_, isMotion, idx);
                        end
                    end
                    return;
                end
                if this.lastNear ~=0 && ~isempty(this.mouseEar) && runFnc
                    feval(this.mouseEar, ax_, isMotion, 0);
                    this.lastNear=0;
                else
                    if DoubleClicker.isMultiSelect
                        this.app.closeToolTip;
                    end
                end
                if this.DEBUG1
                    fprintf('%d/%d (nothing)\n', rx, ry);
                end
            catch
            end
            idx=0;
        end
        
        function [fig_, ax_]=newFig(this, name, standardBars, ...
                visibleOrLocate, property, properties, delayComboBox)
            if nargin<7
                delayComboBox=false;
                if nargin<6
                    properties=[];
                    if nargin<5
                        property=[];
                        if nargin<4
                            visibleOrLocate=true;
                            if nargin<3
                                standardBars=false;
                                if nargin<2
                                    name=sprintf( ...
                                        '2D composite view (of %d dimensions)', ...
                                        this.nDimensions);
                                end
                            end
                        end
                    end
                end
            end
            if isempty(properties)
                properties=BasicMap.Global;
            end
            [fig_, this.tb, personalized, ~]=Gui.Figure( ...
                ~standardBars, property, properties, [], ~standardBars);
            set(fig_, 'UserData', this, 'Name', name);
            if ~standardBars
                if ~personalized
                    P=get(fig_, 'position');
                    set(fig_, 'position', [P(1) P(2) P(3)*1.2 P(4)*1.1]);
                end
                this.btnWrench=ToolBarMethods.addButton( ...
                    this.tb, 'wrench.png', ...
                    'See options', @(h,e)wrench(this, h, e));
                this.btnWrench.setVisible(false);
                n2Tip='See median marker expressions';
                this.btnHeatMap=...
                    ToolBarMethods.addButton(this.tb,'heatMapHot.png', n2Tip, ...
                    @(h,e)n2Msg);
                if ~isempty(this.fncToolBar)
                    feval(this.fncToolBar, this);
                end
                drawnow;
                if ~delayComboBox
                    this.setComboBoxes(fig_);
                end
            end
            ax_=Gui.GetOrCreateAxes(fig_);
            set(ax_, 'Yticklabel', []);
            set(ax_, 'Xticklabel', []);
            this.fig=fig_;
            this.ax=ax_;
            this.setMouse;
            set(this.fig, 'keyPressFcn', ...
                @(src, evd)keyDown(this, src, evd));
            if islogical(visibleOrLocate) && visibleOrLocate
                Gui.SetFigVisible(fig_);
            elseif iscell(visibleOrLocate)
                SuhWindow.Follow(fig_, visibleOrLocate);
                SuhWindow.SetFigVisible(fig_)
            end
            function n2Msg
                this.showMarkerHeatMap;
            end
        end
        
        function svg(this)
            if ~isempty(this.fncSvg)
                feval(this.fncSvg);
            else
                msg('Not yet implemented...');
            end
        end
        
        
        function reconfigureDimensions(this)
            if isempty(this.dimensionName)
                dims='dimensions';
            else
                dims=this.dimensionName;
            end
            if isempty(this.dimensionItemName)
                item='Item';
            else
                item=this.dimensionItemName;
            end
            [filtered, idxs, cancelled]=FilterDlg.Ask(...
                this.allMeasurementNames, this.measurementNames, ...
                'south east+', 'Dimension selector', ['Configure ' dims ' for display'], ...
                item);
            if ~cancelled
                if length(filtered)<2
                    msgBox('At least 2 dimensions are needed');
                    return;
                end
                this.redoDistances(filtered, idxs)
                this.refresh;
            end
        end

        function redoDistances(this, filtered, idxs)
            if nargin<3
                idxs=StringArray.IndexesOf(this.allMeasurementNames, filtered);
            end
            this.measurementNames=filtered;
            this.measurementIdxs=idxs;
            this.measurements=this.allMeasurements(:, idxs);
            this.rawMeasurements=this.allRawMeasurements(:,idxs);
            [~, this.nDimensions]=size(this.measurements);
            this.distances=MatBasics.PDist2Self(this.measurements);
            [this.XY, this.eigvals] = cmdscale(this.distances);
            this.getQuality
        end
        
        function quality=getQuality(this)
            quality=[this.eigvals this.eigvals/max(abs(this.eigvals))]';
        end
        
        function wrench(this, h, ~)
            if ~isempty(this.fncRightClick)
                feval(this.fncRightClick, this, 0, -1, -1, [], 0, h);                
            end
        end
        
        function setComboBoxes(this, fig_)
            hStart=['<html>' this.app.smallStart];
            hEnd=[this.app.smallEnd '</html>'];
            this.cb0=ToolBarMethods.addComboBox([], {...
                [hStart 'No names' hEnd], ...
                [hStart 'Use legend' hEnd], ...
                [hStart 'Use boxes' hEnd], ...
                [hStart 'Boxes+arrows' hEnd], ...
                [hStart 'Order legend by F' hEnd]}, [], true, true, 4, ...
                sprintf('Show names of %s ', this.thingName));
            if isempty(this.propSeeText)
                prop='ndSeeText';
            else
                prop=this.propSeeText;
            end
            this.textPolicy=this.app.getNumeric(prop, 1);
            idx_=this.app.getNumeric(prop, this.textPolicy);
            if idx_==2
                if this.multiProps.is('ndSeeArrow', false)
                    idx_=3;
                end
            end
            this.cb0.setSelectedIndex(idx_);
            set(handle(this.cb0, 'CallbackProperties'), 'ActionPerformedCallback',...
                @(h,e)cb(this, h,e, 0));
            if ~isempty(fig_)
                jp=javaObjectEDT('javax.swing.JPanel');
                if ~this.useComboBoxes
                    jl=javaObjectEDT('javax.swing.JLabel', ...
                        ['<html>' this.app.smallStart '<b>&nbsp;&nbsp;Names:</b>' this.app.smallEnd '</html>']);
                    jp.add(jl);
                end
                jp.add(this.cb0);
                ToolBarMethods.addComponent(this.tb, jp);
                if this.useComboBoxes
                    if size(this.allMeasurements,2)>2
                        try
                            ToolBarMethods.addButton(this.tb, ...
                                'smallGenie.png', 'See dropdowns', ...
                                @(h,e)tinker(this));
                            CytoGate.Get;
                            ToolBarMethods.addButton(this.tb, 'table.gif', ...
                                'Configure dimensions', ...
                                @(h,e)reconfigureDimensions(this));
                        catch
                        end
                    end
                end
            end 
           Gui.AddSvgToToolBar(this.fig, this.tb);
        end
        
        function cb(this, h, e, dim)
            disp(h);
            disp(e);
            disp(dim);
            idx=h.getSelectedIndex;
            if isempty(this.propSeeText)
                prop='ndSeeText';
            else
                prop=this.propSeeText;
            end
            
            switch dim
                case 1
                    this.setSymbolDimension(idx);
                case 2
                    this.setFaceColorDimension(idx);
                case 3
                    this.setEdgeColorDimension(idx);
                case 4
                    this.setTextFaceColorDimension(idx);
                    this.cb0.setSelectedIndex(2);
                    this.app.set(prop, '2');
                case 5
                    this.setTextEdgeColorDimension(idx);
                    this.cb0.setSelectedIndex(2);
                    this.app.set(prop, '2');
                otherwise
                    disp('Re-configuring naming');
                    if idx==3
                        this.textPolicy=2;
                        this.doArrows=true;
                        this.multiProps.set('ndSeeArrow', 'true')
                    else
                        if idx==2
                            this.doArrows=false;
                            this.multiProps.set('ndSeeArrow', 'false')
                        end
                        this.textPolicy=idx;
                    end
                    this.app.set(prop, num2str(this.textPolicy));
            end
            this.refresh;
        end
        function tinker(this)
            ss=this.measurementNames;
            fmt=['<html>' this.app.smallStart '<b>%s %s</b>' this.app.smallEnd '</html>'];
            s=[sprintf(fmt, this.thingName, 'color') ss];
            cb1=ToolBarMethods.addComboBox([], s, ...
                @(h,e)cb(this,h,e, 2), true, true, 1, ...
                sprintf('Alter the color for %s symbols', this.thingName));
            s=[sprintf(fmt, this.thingName, 'border') ss];
            cb2=ToolBarMethods.addComboBox([], s, ...
                @(h,e)cb(this,h,e, 3), true, true, 1, ...
                sprintf('Alter the border color for %s symbols', this.thingName));
            s=[['<html>' this.app.smallStart 'Shape' this.app.smallEnd '</html>'] ss];
            cb3=ToolBarMethods.addComboBox([], s, ...
                @(h,e)cb(this,h,e, 1), true, true, 1, ...
                sprintf('Alter the shapes of %s ', this.thingName));
            jl=javaObjectEDT('javax.swing.JLabel', ...
                ['<html>' this.app.smallStart '<b>&nbsp;&nbsp;Configure:</b>' this.app.smallEnd '</html>']);
            s=[sprintf(fmt, 'text', 'color') ss];
            cb4=ToolBarMethods.addComboBox([], s, ...
                @(h,e)cb(this,h,e, 4), true, true,1, ...
                sprintf('Alter the background color for %s text boxes', this.thingName));
            
            s=[sprintf(fmt, 'text', 'border') ss];
            cb5=ToolBarMethods.addComboBox([], s, ...
                @(h,e)cb(this, h,e, 5), true, true, 1, ...
                sprintf('Alter the border color for %s text boxes', this.thingName));
            jl=javaObjectEDT('javax.swing.JLabel', ...
                ['<html>' this.app.smallStart '<b>&nbsp;&nbsp;Configure:</b>' this.app.smallEnd '</html>']);
            
            jp=javaObjectEDT('javax.swing.JPanel');
            jp.add(jl);
            jp.add(cb1);
            jp.add(cb2);
            jp.add(cb3);
            jp.add(cb4);
            jp.add(cb5);
            msg(jp,0);
        end
            
        function refresh(this)
            this.enable(false);
            if ~isempty(this.btnGenie)
                set(this.btnGenie, 'visible', 'off');
            end
            busyPnl=[];
            try
            busyPnl=BusyPnl(this.fig, 'Refreshing');
            busyPnl.busyJ.stop
            catch
            end
            pu=PopUp('Updating display', 'center', 'Patience...');
            this.initMultiplePlots;
            if isempty(this.fncRefresh)
                this.see;
            else
                feval(this.fncRefresh);
            end
            if ~isempty(busyPnl)
                busyPnl.stop;
            end
            pu.close;
            if ~isempty(this.btnGenie)
                set(this.btnGenie, 'visible', 'on');
            end
            this.enable(true);
            this.updateMrkrExpr;
        end
    end
    
end