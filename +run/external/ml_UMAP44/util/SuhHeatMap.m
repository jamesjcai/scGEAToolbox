classdef SuhHeatMap < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(Constant)
        PROP_COLOR_SCHEME='heatMapScheme';
        COLOR_SCHEMES={...
            'parula', 'jet', 'hot', 'cool', 'spring', 'summer', ...
            'autumn', 'winter', 'gray'};
        COLOR_SCHEME=3;
        ORDER= {...
            'Distance from north east',...
            'Distance from north',...
            'Distance from north west',...
            'Distance from west',...
            'Distance from center',...
            'Distance from east',...
            'Distance from south east',...
            'Distance from south',...
            'Distance from south west',...
            'Ascending', ...
            'Descending', ...
            'Highest max',...
            'Highest average',...
            'Lowest max',...
            'Lowest average'}
    end
    
    methods(Static)
        function p=DefineArgs()
            p = inputParser;
            addParameter(p, 'rawMeasurements', [],  @isMeasurement);
            addParameter(p, 'measurements', [],  @isMeasurement);
            addParameter(p, 'measurementNames', gcf, @(x)isNames(x));
            addParameter(p, 'freqs', [],  @isMeasurement);
            addParameter(p, 'names', {}, @isNames);
            addParameter(p, 'parentFig', gcf, @(x)Gui.IsFigure(x));
            addParameter(p, 'subsetSymbol', [], ...
                @(x)isCallback(x) || Args.IsStrings(x) ...
                || isnumeric(x));
            addParameter(p, 'rowClickedCallback', [], @(x)isCallback(x));
            addParameter(p, 'columnClickedCallback',[], @(x)isCallback(x));
            addParameter(p, 'cellClickedCallback', [], @(x)isCallback(x));
            addParameter(p, 'rowHeaderEnteredCallback', [], @(x)isCallback(x));
            addParameter(p, 'ax', [], @ishandle);
            addParameter(p, 'ignoreScatter', true, @islogical);
            addParameter(p, 'closeFnc', [], @isCallback);
            addParameter(p, 'windowTitleSuffix', ' ', @ischar);
            addParameter(p, 'subsetName', 'Subset', @ischar);
            addParameter(p, 'rowClickAdvice', '', @ischar);
            addParameter(p, 'columnClickAdvice', '', @ischar);
            addParameter(p, 'cellClickAdvice', '', @ischar);
            
            function ok=isCallback(x)
                ok=isequal('function_handle', class(x));
            end
            
            function ok=isNames(x)
                ok=Args.IsStrings(x);
            end
            
            function ok=isMeasurement(x)
                ok=false;
                if isnumeric(x)
                    ok=true;
                end
            end
        end
        
         function [args, errors, argued, unmatchedArgs, ...
                argsObj]=GetArgs(varargin)
            try
            [args, argued, unmatchedArgs, argsObj]=...
                Args.NewKeepUnmatched(...
                SuhHeatMap.DefineArgs(), varargin{:});
            errors=false;
            nFreqs=length(args.freqs);
            nNames=length(args.names);
            if nFreqs ~= nNames
                error('# of names(%d) ~= freqs(%d)', ...
                    nFreqs, nNames);
            end
            nRaw=length(args.rawMeasurements);
            nMea=length(args.measurements);
            if nMea~=nRaw 
                error('# of measurements raw(%d) ~= display(%d)',...
                    nMea, nRaw);
            end
            catch ex
                Gui.MsgException(ex, 'Bad arguments', 9, BasicMap.Global);
            end
         end
        
        function [top, cbColor, cbDist]=DropDowns(fncColor, propDist, ...
                fncDist, dfltDist, prependDist, whereColor, whereDist, ...
                sfx, useDistanceOrder)
            if nargin<9
                useDistanceOrder=true;
                if nargin<8
                    sfx='';
                    if nargin<7
                        whereDist='North';
                        if nargin<6
                            whereColor='South';
                        end
                    end
                end
            end
            cbColor=SuhHeatMap.ColorSchemeDropDown(fncColor);
            cbDist=SuhHeatMap.DistanceDropDown(propDist, ...
                fncDist, dfltDist, prependDist, useDistanceOrder);
            if useDistanceOrder
                top=Gui.BorderPanel;
                topWest=Gui.FlowLeftPanel(4,0, ...
                    Html.WrapSmallBold('Sort: '), cbDist);
                top.add(topWest, whereDist);
                topEast=Gui.FlowLeftPanel(0,0,...
                    Html.WrapSmallBold([sfx ': ']), cbColor);
                top.add(topEast, whereColor);
            else
                top=Gui.FlowLeftPanel(4,0,...
                Html.WrapSmallBold('Sort: '), cbDist,...
                Html.WrapSmallBold([sfx ': ']), cbColor);
            end
        end
        
        function cb=ColorSchemeDropDown(fnc)
            app=BasicMap.Global;
            prop=SuhHeatMap.PROP_COLOR_SCHEME;
            cs=SuhHeatMap.COLOR_SCHEMES;
            N=length(cs);
            lbls=cell(1,N);
            for i=1:N
                lbls{i}=['<html>' ...
                    Html.ImgXy([cs{i} '.png'], [], .33) ...
                    '</html>'];
            end
            cb=Gui.Combo(lbls, SuhHeatMap.COLOR_SCHEME, prop, [], fnc);
            cb.setMaximumRowCount(N);
        end
        
        function cb=DistanceDropDown(prop, fnc, dflt, prepend, useDistance)
            if nargin<4
                prepend={};
            end
            if useDistance
                opts=Html.WrapSmallBoldCell(...
                    [prepend SuhHeatMap.ORDER]);
            else
                opts=Html.WrapSmallBoldCell(...
                    [prepend SuhHeatMap.ORDER(end-5:end)]);
            end
            cb=Gui.Combo(opts, dflt, prop, [], fnc);
            cb.setMaximumRowCount(length(opts));
        end
        
        function [I, low2High]=SortDists(ax_, xy, distFromWhere, xl, yl)
            if nargin<3
                distFromWhere='east';
            end
            distFromWhere=char(edu.stanford.facs.swing.Basics.RemoveXml(...
                distFromWhere));
            if String.StartsWith(distFromWhere, 'Distance from ')
                distFromWhere=distFromWhere(15:end);
            end
            if nargin<5
                yl=ylim(ax_);
            end            
            if strcmpi('reverse', get(ax_, 'YDir'))
                yl=[yl(2) yl(1)];
            end
            if nargin<4
                xl=xlim(ax_);
            end            
            if strcmpi('reverse', get(ax_, 'XDir'))
                xl=[xl(2) xl(1)];
            end
            midX=(xl(1)+(xl(2)-xl(1)/2));
            midY=(yl(1)+(yl(2)-yl(1)/2));
            if strcmpi('north east',  distFromWhere)
                pt=[xl(2) yl(2)];
            elseif strcmpi('north',  distFromWhere)
                pt=[midX yl(2)];
            elseif strcmpi('north west',  distFromWhere)
                pt=[xl(1) yl(2)];
            elseif strcmpi('south west',  distFromWhere)
                pt=[xl(1) yl(1)];
            elseif strcmpi('south',  distFromWhere)
                pt=[midX yl(1)];
            elseif strcmpi('south east',  distFromWhere)
                pt=[xl(2) yl(1)];
            elseif strcmpi('east',  distFromWhere)
                pt=[xl(2) midY];
            elseif strcmpi('west',  distFromWhere)
                pt=[xl(1), midY];
            elseif strcmpi('center',  distFromWhere)
                pt=[midX, midY];
            end
            D=pdist2(pt, xy);
            [low2High,I]=sort(D);
        end
        
        function [jd, html]=Msg(args, showNoMatch, rowIdxs, ...
                colIdxs, where, ttl)
            if nargin<6
                ttl='Median marker expressions';
                if nargin<5
                    where='east+';
                end
            end
            if nargin==1
                [html, tp]=SuhHeatMap.Html(args);
            elseif nargin==2
                [html, tp]=SuhHeatMap.Html(args, [], showNoMatch);
            elseif nargin==3
                [html, tp]=SuhHeatMap.Html(args, [], showNoMatch, rowIdxs);
            elseif nargin==4
                [html, tp]=SuhHeatMap.Html(args, [], showNoMatch, rowIdxs, colIdxs);
            end
            msg(tp,0,where,ttl, 'none');
        end
        
        function [bp, tp]=TextPane(args, ax_, rowHdrs)
            app=BasicMap.Global;
            tp=Gui.TextPane;
            bp=Gui.BorderPanel(1,5);
            bp.add(Gui.Scroll(tp), 'Center');
            lbl=Gui.Label(' ');
            lbl.setHorizontalAlignment(javax.swing.JLabel.CENTER);
            bp.add(lbl, 'South');
            hj= handle(tp,'CallbackProperties');
            hovering=-1;
            tipIdx=0;
            set(hj,'HyperlinkUpdateCallback', ...
                @(h,e)LinkCallback(h, e, args, ax_, rowHdrs,lbl));
            isInt=range(args.rawMeasurements)>10;
            
            function LinkCallback(h, e, args, ax_, rowHdrs, jl)
                
                description = char(e.getDescription); % URL stri
                et=char(e.getEventType);
                idxs=str2num(description(2:end)); %#ok<ST2NM> 
                mouseEvent=e.getInputEvent;
                
                switch char(et)
                    case char(e.getEventType.EXITED)
                        jl.setText(' ');
                        if ~isempty(ax_)
                            BasicMap.Global.closeToolTip;
                        end
                    case char(e.getEventType.ACTIVATED)
                        if idxs(1)>0
                            rowName=args.names{idxs(1)};
                        else
                            rowName='column header';
                        end
                        if idxs(3)>0
                            colName=args.measurementNames{idxs(3)};
                        else
                            colName='row header';
                        end
                        fprintf(['Clicked on "' description ...
                            '" row="%s"' ' column="%s"\n'], ...
                            rowName, colName);
                        if idxs(1)>0
                            if idxs(3)==0
                                if ~isempty(args.rowClickedCallback)
                                    feval(args.rowClickedCallback, ...
                                        ax_, idxs(1));
                                end
                            elseif ~isempty(args.cellClickedCallback)
                                feval(args.cellClickedCallback, ...
                                    ax_, idxs(1), idxs(3));
                            end
                        elseif ~isempty(args.columnClickedCallback)
                            feval(args.columnClickedCallback, ...
                                ax_, idxs(3));
                        end
                    case char(e.getEventType.ENTERED)
                        if length(idxs)~=4
                            if isempty(ax_)
                                jl.setText(' ');
                            end
                            return;
                        end
                        if idxs(3)>0
                            if idxs(1)>0
                                if ~isInt
                                    numb=String.encodeRounded(...
                                        args.rawMeasurements(idxs(1), idxs(3)),2);
                                else
                                    numb=String.encodeInteger(...
                                        args.rawMeasurements(idxs(1), idxs(3)));
                                end
                                tip=sprintf('%s=<b>%s</b>', ...
                                    String.ToHtml(strrep(...
                                    args.measurementNames{idxs(3)}, '\bf','')), ...
                                    numb);
                                
                            else
                                tip=String.ToHtml(strrep(...
                                    args.measurementNames{idxs(3)}, '\bf',''));
                                if isfield(args, 'columnClickAdvice')...
                                        && ~isempty(args.columnClickAdvice)
                                    if contains(args.columnClickAdvice,'%s')
                                        %is sprintf string
                                        tip=sprintf(...
                                            args.columnClickAdvice, ...
                                            ['<font color="blue"><b>'...
                                            tip '</b></font>']);
                                    else
                                        tip=[tip ' ' args.columnClickAdvice];
                                    end
                                end
                            end
                        else
                            tip='';
                        end
                        if idxs(1)>0
                            if idxs(3)>0
                                tip=[tip ' for subset '];
                            end
                            tip=sprintf('%s "%s"', ...
                                tip, rowHdrs{idxs(1)});
                            if idxs(3)>0
                                if isfield(args, 'cellClickAdvice')...
                                        && ~isempty(args.cellClickAdvice)
                                    tip=[tip ' ' args.cellClickAdvice];
                                end
                            else
                                if isfield(args, 'rowClickAdvice')...
                                        && ~isempty(args.rowClickAdvice)
                                    tip=[tip ' ' args.rowClickAdvice];
                                end
                            end
                        else
                        end
                        %disp('link hover enter');
                        jl.setText(['<html><b>' args.app.smallStart  ...
                            tip args.app.smallEnd '</b></html>']);
                        
                        if idxs(3)==0 &&...
                                ~isempty(args.rowHeaderEnteredCallback) ...
                                && idxs(1)~=tipIdx
                            tipIdx=idxs(1);
                            if idxs(1)>0 && ~isempty(ax_)
                                if hovering ~= idxs(1)
                                    hovering=idxs(1);
                                    [tip, botPnl]=feval(...
                                        args.rowHeaderEnteredCallback, ...
                                        ax_, idxs(1));
                                    xOffset=mouseEvent.getX+4;
                                    yOffset=mouseEvent.getY-11;
                                    app.showToolTip(h, tip, xOffset, ...
                                        yOffset, 5, botPnl);
                                end
                            end
                        end
                end
            end
        end
        
        function schemeName=DefaultColorScheme
            idx=BasicMap.Global.getNumeric(...
                SuhHeatMap.PROP_COLOR_SCHEME,...
                SuhHeatMap.COLOR_SCHEME);
            schemeName=SuhHeatMap.COLOR_SCHEMES{idx};
        end
        
        function [map, schemeName]=ColorScheme
            schemeName=SuhHeatMap.DefaultColorScheme;
            map=feval(schemeName);
        end
        
         function [html, bp, tp]=Html(args, ax_, showNoMatch, rowIdxs,...
                colIdxs, ignoreScatter, forBrowser, useColNums, ...
                scaleName, border, padding, spacing, rowLabel, ...
                columnLabel, imgFldr)
            if nargin<15
                imgFldr=[];
                if nargin<14
                    columnLabel=[];
                    if nargin<13
                        rowLabel=[];
                        if nargin<12
                            spacing='0';
                            if nargin<11
                                padding='0';
                                if nargin<10
                                    border='0';
                                    if nargin<9
                                        scaleName='Low to high expression';
                                        if nargin<8
                                            useColNums=true;
                                            if nargin<7
                                                forBrowser=false;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            [J, scheme]=SuhHeatMap.ColorScheme;
            nJ=size(J,1);
            [R,C]=size(args.measurements);
            names=args.names;
            app=args.app;
            N=length(names);
            hasHtml=String.StartsWith(names{1}, '<html>');
            for i=1:N
                if any(strfind(names{i}, '\bf'))
                    names{i}=app.getMatchName(String.ToHtml(names{i}), ...
                        false);
                elseif hasHtml
                    names{i}=Html.remove(names{i});
                else
                    names{i}=String.ToHtml(names{i});
                end
            end
            if nargin<6
                ignoreScatter=false;
                if nargin<5
                    colIdxs=[];
                    if nargin<4
                        rowIdxs=[];
                        if nargin<3
                            showNoMatch=true;
                            if nargin<2
                                ax_=[];
                                forBrowser=true;
                            end
                        end
                    end
                else
                    C=length(colIdxs);
                end
            end
            mNames=args.measurementNames;
            N=length(mNames);
            for i=1:N
                if String.Contains(mNames{i}, '\bf')
                    mNames{i}=[strrep(String.ToHtml(mNames{i}), '\bf', '<b>') '</b>'];
                else
                    mNames{i}=String.ToHtml(mNames{i});
                end
            end
            if isempty(colIdxs)
                colIdxs=1:C;
            end
            dataColIdxs=SuhHeatMap.FilterCols(...
                args, ignoreScatter, colIdxs);
            maxC=C;
            C=length(dataColIdxs);
            if isempty(rowIdxs)
                rowIdxs=1:R;
            else
                R=length(rowIdxs);
            end
            if forBrowser
                sm1='<small>';
                sm2='</small>';
            else
                sm1=args.app.smallStart;
                sm2=args.app.smallEnd;
            end
            if forBrowser
                tmpBorder='0';
            else
                tmpBorder=border;
            end
            if isempty(rowLabel) && isempty(columnLabel)
                if isfield(args, 'subsetName')
                    sn=args.subsetName;
                else
                    sn='Subset';
                end
                rowLabel=[sm1 '&nbsp;&nbsp;<font color="blue">' ...
                    Html.Img('downArrow.png') sn ...
                    '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'...
                    '&nbsp;&nbsp;&nbsp;<i>'...
                    'Measurements' Html.Img('rightArrow.png')...
                    '</i></font>' sm2 '<hr>'];
                columnLabel='';
            elseif ~isempty(columnLabel)
                columnLabel=['<tr><th></th><th colspan="' num2str(C+2)...
                    '">' columnLabel '<hr></th></tr>'];
            end
            sb=java.lang.StringBuilder(Html.RotatedStyle);
            sb.append(['<STYLE>A {text-decoration: none;} </STYLE>'...
                '<STYLE>br {line-height: .5;}</STYLE>'...
                '<table border="']);
            sb.append(tmpBorder);
            sb.append('" cellpadding="');
            sb.append( padding);
            sb.append( '" cellspacing="');
            sb.append( spacing );
            sb.append('" bgcolor="white"><thead>');
            sb.append(columnLabel );
            sb.append('<tr><th align="left" valign="bottom">');
            sb.append(rowLabel);
            sb.append('</th>');
            if ~forBrowser
                sb.append('<th valign="bottom"> ');
                sb.append(Html.Vertical('Freq %', sm1, sm2));
                sb.append('<hr></th><th valign="bottom"><hr></th>');
            else
                sb.append(Html.Rotate('Frequency %'));
                sb.append('<th></th>');
            end
            for col=1:C
                dataColIdx=dataColIdxs(col);
                tip=['0 0 ' num2str(dataColIdx) ' ' num2str(col)];
                colNum='';
                if ~forBrowser
                    a1=['<a href="#' tip '">' sm1];
                    a2=[sm2 '</a>'];
                    if useColNums
                        colNum=[a1 '<i>' num2str(dataColIdx) ...
                            '</i>' a2 '<hr>'];
                    end
                    sb.append('<th valign="bottom"> ');
                    sb.append(Html.Vertical(mNames{dataColIdx}, a1, a2));
                    sb.append(   colNum);
                    sb.append( '</th>');
                else
                    sb.append(Html.Rotate(mNames{dataColIdx}));
                end
            end
            sb.append('</tr></thead>');
            for row=1:R
                dataRowIdx=rowIdxs(row);
                if ~showNoMatch
                    if args.matchIdxs(dataRowIdx)==0
                        continue;
                    end
                end
                tipStart=[num2str(dataRowIdx) ' ' num2str(row) ' '];
                if ismethod(args, 'subsetSymbol')
                    subsetColor=[args.subsetSymbol(dataRowIdx) ...
                        '&nbsp;&nbsp;'];
                elseif iscell(args.subsetSymbol)
                    subsetColor=[args.subsetSymbol{dataRowIdx}...
                        '&nbsp;&nbsp;'];
                elseif isnumeric(args.subsetSymbol)
                    clr=args.subsetSymbol(dataRowIdx, :);
                    fs=floor(5+args.freqs(dataRowIdx)*55);
                    str=['<font  ' Gui.HtmlHexColor(clr) ...
                        '>&bull;</font>'];
                    subsetColor=['<font size="' num2str(fs) ...
                        '">' str '</font>'];                    
                else
                    subsetColor=[feval(args.subsetSymbol, dataRowIdx) ...
                        '&nbsp;&nbsp;'];
                end
                freq=args.freqs(dataRowIdx);
                perc=String.encodePercent(freq, 1, 1);
                perc(end)=' ';
                tip=[tipStart ' 0 0 ' ];
                sb.append('<tr><td><a href="#');
                sb.append(tip);
                sb.append('">');
                sb.append(subsetColor );
                name=names{dataRowIdx};
                sb.append('<b>');
                sb.append(sm1);
                sb.append(name);
                if hasHtml
                    imgIdx=String.IndexOf(name, '<img ');
                    if imgIdx>0
                        sb.append((1:imgIdx-1));
                        sb.append('&nbsp;&nbsp;</a>');
                        sb.append(name(imgIdx:end));
                    else
                        sb.append('&nbsp;&nbsp;</a>');
                    end
                else
                    sb.append('&nbsp;&nbsp;</a>');
                end
                sb.append(sm2);
                sb.append('</b></td><td align="right">');
                sb.append(perc);
                sb.append('</td><td>&nbsp;&nbsp;&nbsp;</td>');
                for col=1:C
                    dataColIdx=dataColIdxs(col);
                    pk=args.measurements(dataRowIdx, dataColIdx);
                    if isnan(pk)
                        sb.append('<td></td>');
                        continue;
                    end
                    jIdx=floor(pk*nJ);
                    try
                        c=J(jIdx,:);
                    catch ex
                        disp(ex);
                        if jIdx<1
                            c=J(1,:);
                        else
                            c=J(end,:);
                        end
                    end
                    tip=[tipStart num2str(dataColIdx) ' ' num2str(col)];
                    sb.append('<td bg');
                    sb.append(Gui.HtmlHexColor(c));
                    sb.append('><a href="#');
                    sb.append(tip);
                    sb.append('">&nbsp;&nbsp;&nbsp;&nbsp;</a></td>');
                end
                sb.append('</tr>');
            end
            schemePng=[scheme '.png'];
            if ~isempty(imgFldr)
                fromPng=fullfile(BasicMap.Global.contentFolder, schemePng);
                [~, imgSubFldr]=fileparts(imgFldr);
                toPng=fullfile(imgFldr, schemePng);
                File.Copy(fromPng, toPng);
                img=Html.ImgXy(schemePng, imgSubFldr, .71, forBrowser) ;
            else
                img=Html.ImgXy(schemePng, [], .71, forBrowser) ;
            end
            sb.append('<tr><td align="center" colspan="');
            sb.append(num2str(maxC+3));
            sb.append('"><br><br><i><b>');
            sb.append(sm1)
            sb.append(scaleName);
            sb.append(sm2);
            sb.append('</b></i><br>' );
            sb.append(img);
            sb.append('</td></tr></table>');
            html=char(sb.toString);
            if nargout>1
                [bp, tp]=SuhHeatMap.TextPane(args, ax_, names);
                tp.setText(['<html>' html '<hr><br><br></html>']);
            end
        end
 
        function I=SortByMeasurement(key, args, rowIdxs)
            [~,C]=size(args.measurements);
            colIdxs=SuhHeatMap.FilterCols(args, args.ignoreScatter, 1:C);
            if nargin<3
                data=args.measurements(:, colIdxs)';
            else
                data=args.measurements(rowIdxs, colIdxs)';
            end
            if length(colIdxs)>1
                if endsWith(key, 'est max')
                    data=max(data);
                else
                    data=mean(data);
                end
            end
            if startsWith(key, 'Highest')
                [~,I]=sort(data, 'descend');
            else
                [~,I]=sort(data, 'ascend');
            end
        end
        function dataColIdxs=FilterCols(args, ignoreScatter, colIdxs)
            dataColIdxs=[];
            if ignoreScatter
                cns=args.measurementNames;
                nCns=length(cns);
                ignoreIdxs=[];
                for j=1:nCns
                    if String.StartsWith(cns{j}, 'SSC-') || ....
                            String.StartsWith(cns{j}, 'FSC-')
                        ignoreIdxs(end+1)=j;
                    end
                end
            else
                ignoreIdxs=[];
            end
            ignoring=~isempty(ignoreIdxs);
            selectedMarkers=[];
            try
                selectedMarkers=args.selectedMarkers;
            catch
            end
            
            C=length(colIdxs);
            for col=1:C
                dataColIdx=colIdxs(col);
                if ignoring
                    if any(dataColIdx==ignoreIdxs)
                        continue;
                    end
                end
                if ~isempty(selectedMarkers)
                    if ~any(dataColIdx==selectedMarkers)
                        continue;
                    end
                end
                dataColIdxs(end+1)=dataColIdx;
            end
        end

        function [jd, tp1]=New(varargin)
            jd=[];
            firstPack=true;
            tp1=[];
            [args, errors]=SuhHeatMap.GetArgs(varargin{:});
            if errors>0
                return;
            end
            args.app=BasicMap.Global;
            if isempty(args.subsetSymbol)
                args.subsetSymbol=@symb;
            end
            [~,C]=size(args.measurements);
            options=cell(1,C);
            for col=1:C
                options{col}=['<html>' args.app.smallStart ...
                    '<font color="blue"><i>' ...
                    num2str(col) '.</i></font> <b>' ...
                    args.measurementNames{col} '</b>'...
                    args.app.smallEnd '</html>'];
            end
            showNoMatch=true;
            heatMap=Gui.BorderPanel(10,1);
            [top, ~, cb]=SuhHeatMap.DropDowns(...
                @(h,e)setColorScheme, 'originalTreeHeatMapOrder', ...
                @(h,e)reorder, 1, {'Original order'}, ...
                'South', 'North', 'Low to high', false);
            refreshBtn=Gui.NewBtn(Html.WrapSm('Refresh'), @(h,e)refresh(),...
                'Click to synchronize HeatMap to marker selections ', ...
                'refresh.png');
            refreshBtn.setEnabled(false);
            heatMap.add(Gui.BorderPanel([],1,1, 'Center', top, ...
                'East', Gui.Panel(refreshBtn)), 'North');
            jd=[];
            bp1=[];
            maxN=0;
            I=[];
            setMap;
            nItems=18; 
            if maxN<7
                nItems=15; % with column height and color legend 15 is fine
            elseif maxN<12
            elseif maxN<20
                nItems=25;
            else
                nItems=30; % 28 colors for BD instruments
            end
            figure(args.parentFig);
            browseBtn=Gui.NewBtn(Html.WrapSm('Browse'), ...
                @(h,e)browse(),...
                'Click to see web page ', ...
                'world_16.png');
            sw=Gui.Panel(browseBtn);
            filteredDataColIdxs=SuhHeatMap.FilterCols( ...
                args, args.ignoreScatter, 1:C);
            options=options(filteredDataColIdxs);
            cdI=edu.stanford.facs.swing.MarkerSorter.Sort(options);
            opts=options(cdI);
            idxsToDo=[];
            [~,~,jd]=mnuMultiDlg(struct(...
                'allMsg', '<html><b>All</b><hr></html>',...
                'noCancel', true, 'closeFnc', args.closeFnc, ...
                'where', 'west++','msg', '', 'icon', 'none', ...
                'modal', false, 'checkBoxFnc',...
                @(h, e, idxs, checkBoxes)markerSelected(idxs)),...
                ['HeatMap' args.windowTitleSuffix], opts, ...
                [], false, true, sw, ...
                'south west buttons', heatMap,'West',  nItems);
            firstPack=SuhHeatMap.Repack(jd, firstPack);
            
            function reorder
                jd.setEnabled(false);
                setMap;
                jd.setEnabled(true);
            end
            
            
            function markerSelected(idxs)
                idxsToDo=idxs;
                if isempty(idxs) 
                    refreshBtn.setEnabled(false);
                    return;
                end
                refreshBtn.setEnabled(true);
                edu.stanford.facs.swing.Basics.Shake(refreshBtn, 3);
                BasicMap.Global.showToolTip(refreshBtn, [], 20, 25);
            end

            function refresh
                jd.setEnabled(false);
                Gui.ShowFacs(jd, 'Refreshing HeatMap');
                try
                    args.selectedMarkers=filteredDataColIdxs(cdI(idxsToDo));
                    setMap;
                catch
                end
                Gui.HideBusy(jd);
                jd.setEnabled(true);
                refreshBtn.setEnabled(false);
            end
            
            function setColorScheme
                jd.setEnabled(false);
                setMap;
                jd.setEnabled(true);
            end
            
            function browse()
                figure(args.parentFig);
                pu=PopUp('Creating webpage html', 'west+');
                htmlHeat=SuhHeatMap.Html(args, [], ...
                    showNoMatch, I, [], args.ignoreScatter, true);
                if ~isempty(args.ax) && ~isempty(args.parentFig)
                    pngPath=[tempname '.png'];
                    [imgFldr, f, e]=fileparts(pngPath);
                    pngFile=[f e];
                    F = getframe(args.ax);
                    Image = frame2im(F);
                    imwrite(Image, pngPath);
                    saveas(args.parentFig, pngPath);
                    htmlHeat=['<table><tr><td colspan="2"><br><hr>' ...
                        Html.ImgForBrowser(pngFile, imgFldr)...
                        '</td></tr></table><table cellpadding="14">'...
                        '<tr><td valign="top">' ...
                        htmlHeat '</td></tr></table>'];
                
                end
                Html.Browse(['<html>' htmlHeat '</html>']);
                pu.close;
            end
            
            function str=symb(idx)
                clr=Gui.HslColor(idx, length(args.names));
                fontSize=floor(5+args.freqs(idx)*55);
                str=['<font  ' Gui.HtmlHexColor(clr)...
                    '>&bull;</font>'];
                str=['<font size="' num2str(fontSize) '">' str '</font>'];
            end
            
            function setMap
                subsetOrder=char(cb.getSelectedItem);
                if ~isempty(bp1)
                    heatMap.remove(bp1);
                end
                I=sortXyAndAx(subsetOrder);
                maxN=length(I);
                [~, bp1, tp1]=SuhHeatMap.Html(args, [], ...
                    showNoMatch, I, [], args.ignoreScatter, false);
                heatMap.add(bp1, 'Center');
                firstPack=SuhHeatMap.Repack(jd, firstPack);
            end
            
            function I=sortXyAndAx(distFromWhere) 
                distFromWhere=char(...
                    edu.stanford.facs.swing.Basics.RemoveXml(...
                    distFromWhere));
                N1=length(args.names);
                if startsWith(distFromWhere, 'Ascending')
                    I=sort_(args.names);
                elseif startsWith(distFromWhere, 'Descending')
                    I=sort_(args.names);
                    I=flip(I);
                elseif startsWith(distFromWhere, 'Highest') || ...
                        startsWith(distFromWhere, 'Lowest')
                    I=SuhHeatMap.SortByMeasurement(distFromWhere, args);
                else
                    I=1:N1;
                end
                function I=sort_(strs)
                    out=StringArray.ToLower(strs);
                    N=length(out);
                    for ii=1:N
                        out{ii}=strrep(out{ii}, '\bf', '');
                    end
                    [~, I]=sort(out);
                end
            end
            
        end
        
        function firstPack=Repack(jd, firstPack)
            if ~isempty(jd)
                if ~firstPack
                    dim=jd.getSize;
                end
                jd.pack;
                if ~firstPack
                    dim2=jd.getSize;
                    dim2.height=dim.height;
                    jd.setSize(dim2);
                else
                    firstPack=false;
                end
            end
        end
            
    end
end
