classdef SuhScaler <SuhAbstractClass
%   Class framework to support heterogenous methods for normalizing data
%   scales to 0-1. Helps for flow cytometry where lower small regions of
%   dull signal have significant information lost to clustering methods
%   unless stretched.
%   Inspired by Wayne Moore's Transform GitHub project in Java which lacks
%   fast vectorization.

%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(Constant)
        DEFAULT_DECADES=4.5;%Number of decades in default log or log like scale
        LN_10=log(10); %Constant needed by log and log like scales
        LOG='log';
        LINEAR='linear';
        MILTENYI='miltenyi';
        HYPERLOG='hyperlog';
        ARCSINH='fasinh';
        BIEX='biex'
        LOGICLE='logicle';
    end

    properties(SetAccess=protected)
        T; %Top of scale for all methods
        A; %Bottom of scale
        M; %Decades for methods which need this
        W; %Width for methods needing 2 types of scale (e.e. linear and log)
        J; %A Java transform object for operations not needing MATLAB
           %translation for fast vectorization (e.g. label suggestions)
        bins;
        expLabels=true;
        labels=[];
        positions;
        scaledPositions;
        scaledTicks;
        tickSizes;
        isLog10Scale;
        minDisplay;
        maxDisplay;
        zeroPosition;
        % For scales with negative regions Stephen likes to see negative
        % labels
        seeNegativeLabels=true;
        lastDisplayNeg;
    end
    
    methods
        function this=SuhScaler(T, bins)
            if nargin<2
                bins=0;
            end
            this=this@SuhAbstractClass();
            this.bins=bins;
            this.T=T;
        end
        
        function scaleValues=scale(this, dataValues)
            %Computes a 0-1 scale position for 1 or more data values
            %
            %Input arguments
            %value          vector of values from external scale
            %
            %Output arguments
            %position of value on 0-1 scale
            
            if isempty(this.J)
                this.warnNotImplemented;
                scaleValues=nan;
            else
                if isempty(dataValues)
                    scaleValues=[];
                else
                    scaleValues=this.J.scale(dataValues);
                end
            end
        end
        
        function dataValues=inverse(this, scaleValues)
            %Computes the inverse data value for 1 or more
            %   values on the 0-1 normalized scale
            %
            %Input arguments
            %value    values on this method's normalized scale
            %
            %Output arguments
            %position of value on 0-1 scale
            
            if isempty(this.J)
                this.warnNotImplemented;
                dataValues=nan;
            else
                dataValues=this.J.inverse(scaleValues);
            end
        end
        
        function setSeeNegativeValues(this, ok)
            this.seeNegativeLabels=ok;
            this.labels=[];%cause recompute of labels
        end
        
        %A suggested list of data values suitable for use as 
        %   labels on an axis for this transform. Each data value 
        %   may be turned into a string to be located at the scale 
        %   position of the data.
        function [positions, scaledPositions, labels, zeroPosition, ...
                minDisplay, maxDisplay, scaledTicks, tickSizes, ...
                isLog10Scale]=getDisplayInfo(this)
            if isempty(this.labels)                
                [minDisplay, maxDisplay]=this.getDisplayLimits;
                if isempty(this.zeroPosition)
                    % may be initialized to NaN for scalers like LOG
                    zeroPosition=this.scale(0);
                else
                    zeroPosition=this.zeroPosition;
                end
                adjustForNegative=true;
                if isempty(this.J)
                    this.warnNotImplemented;
                    positions=[];
                    scaledPositions=[];
                    labels={};
                    scaledTicks=[];
                    tickSizes=[];
                    return;
                else
                    try
                        positions=this.J.axisLabels;
                    catch
                        adjustForNegative=false;
                        positions=MatBasics.GetPositions(this.A, this.T);
                    end
                end
                nPositions=length(positions);
                this.lastDisplayNeg=minDisplay;
                mn2=0;
                if this.seeNegativeLabels && adjustForNegative
                    if ~isnan(zeroPosition) ...
                            && minDisplay<zeroPosition
                        top=floor(log10(this.T));
                        mn2=this.inverse(minDisplay);
                        if top<6
                            if mn2<-200
                                mn2=-1000;
                            else
                                mn2=-100;
                            end
                        else
                            mn2=-10000;
                        end
                        mm=log10(abs(mn2))-1;
                        fp=floor(positions);
                        while mm>0
                            if ~isempty(find(fp==10^mm, 1))
                                num=0-10^mm;
                                if isempty(find(fp==num,1))
                                    mn2=[mn2;num];
                                end
                            end
                            mm=mm-1;
                        end
                        scaledMn2=this.scale(mn2(1));
                        if scaledMn2<this.lastDisplayNeg
                            this.lastDisplayNeg=scaledMn2;
                        end
                        if mn2(1)<positions(1)
                            positions=[mn2;positions];
                        else
                            %disp('huh');
                        end
                    end
                end
                N=length(positions);
                labels=cell(1,N);
                if isnan(positions)
                    labels{1}='N/A';
                else
                    if this.expLabels
                        for i=1:N
                            if (positions(i)==1 || positions(i)==-1 ) ...
                                    && nPositions>3
                                labels{i}='';
                            else
                                labels{i}=String.encodePowTen(positions(i));
                            end
                        end
                    else
                        top=floor(log10(this.T));
                        if top<5
                            unit=1;
                            suffix='';
                        elseif top<6
                            unit=1000;
                            suffix='K';
                        else
                            unit=1000000;
                            suffix='M';
                        end
                        for i=1:N
                            p=positions(i);
                            if p>=unit
                                labels{i}=[String.encodeInteger(...
                                    p/unit) suffix];
                            else
                                labels{i}=String.encodeInteger(p);
                            end
                        end
                    end
                end
                last=log10(positions(end));
                isLog10Scale=isequal(last, floor(last));
                if isLog10Scale
                    increments=log10([20 30 40 50 60 70 80 90])-1;
                    increments(end+1)=1;
                    if mn2(1)<0
                        start=0-ceil(log10(abs(mn2(1))));
                    else
                        start=0;
                    end
                    nMajors=(last-start)+1;
                    nIncrements=length(increments);
                    nTicks=nMajors*nIncrements;
                    scaledTicks=zeros(1,nTicks);
                    tickSizes=zeros(1,nTicks);
                    idx=1;
                    for i=start:last
                        for j=1:nIncrements
                            if i<0
                                next=(abs(i)-1)+increments(nIncrements-j+1);
                                num=0-(10^next);
                                if j==1
                                    tickSizes(idx)=2;
                                else
                                    tickSizes(idx)=1;
                                end
                            else
                                next=i+increments(j);
                                if next==0
                                    num=0;
                                else
                                    num=10^next;
                                end
                                if j==nIncrements
                                    tickSizes(idx)=2;
                                else
                                    tickSizes(idx)=1;
                                end
                            end
                            scaledTicks(idx)=this.scale(num);
                            idx=idx+1;
                        end
                    end
                    if ~isempty(find(positions==-1,1))
                        scaledTicks(end+1)=this.scale(-1);
                        tickSizes(end+1)=2;
                    end
                    if ~isempty(find(positions==1,1))
                        scaledTicks(end+1)=this.scale(1);
                        tickSizes(end+1)=2;
                    end

                else
                    increment=(positions(end)-positions(end-1))/5;
                    nTicks=length(positions)*5;
                    scaledTicks=zeros(1,nTicks);
                    tickSizes=zeros(1,nTicks);
                    next=increment;
                    for i=1:nTicks
                        if mod(i, 5)==0
                            tickSizes(i)=2;
                        else
                            tickSizes(i)=1;
                        end
                        scaledTicks(i)=this.scale(next);
                        next=next+increment;
                    end
                end
                this.positions=positions;
                this.scaledPositions=this.scale(positions');
                this.labels=labels;
                this.zeroPosition=zeroPosition;
                this.minDisplay=minDisplay;
                this.maxDisplay=maxDisplay;
                this.scaledTicks=scaledTicks;
                this.tickSizes=tickSizes;
                this.isLog10Scale=isLog10Scale;
            end
            positions=this.positions;
            scaledPositions=this.scaledPositions;
            labels=this.labels;
            zeroPosition=this.zeroPosition;
            minDisplay=this.minDisplay;
            maxDisplay=this.maxDisplay;
            scaledTicks=this.scaledTicks;
            tickSizes=this.tickSizes;
            isLog10Scale=this.isLog10Scale;
        end
        
        function [mn, mx]=getDisplayLimits(~)
            mn=0;
            mx=1;
        end
        
        function inverted=inverseFromBins(this, bin, maxBin)
            inverted=this.inverse(bin/maxBin);
        end

        function bin=unscaledToBins(this, inverted, maxBin)
            bin=SuhScaler.ScaledToBins(this.scale(inverted), maxBin);
        end

        function value=percent(this, percent, dim)
            if isempty(this.A)
                mn=0;
            else
                mn=this.A;
                if mn~=0
                    value=this.A+(percent*(this.T-this.A));
                    if nargin>2
                        prefix=['of "' dim '" ' ];
                    else
                        prefix='';
                    end
                    sp=String.encodePercent(percent);
                    sv=num2str(value);
                    fprintf( ...
                        ['NOTE:  %s %sis %s; forumula: ' ...
                        'A+(%s*(T-A)) where A=%s, T=%s!\n'], ...
                        sp, prefix, sv, sp, num2str(this.A), ...
                        num2str(this.T));
                    return;
                end
            end
            value=percent*(this.T-mn);
        end
    end
    
    methods(Static)
        function bin=ScaledToBins(scaled, maxBin)
            bin=scaled*maxBin;
        end

        function yes=HasJava
            persistent hasJava;
            if isempty(hasJava)
                hasJava=SuhScaler.InitJava;
            end
            yes=hasJava;
        end
       
        function [ok, test]=InitJava
            ok=false;
            try
                test=edu.stanford.facs.transform.Hyperlog(10000, 2.5, 5.0);
            catch
                jar=fullfile(fileparts(mfilename('fullpath')), ...
                    'transform.jar');
                javaaddpath(jar);
            end
            try
                test=edu.stanford.facs.transform.Hyperlog(10000, 2.5, 5.0);
                ok=true;
            catch ex
                ex.getReport
            end
        end
       
        function matrix=Inverse(scalers, matrix)
            [~,C]=size(matrix);
            for c=1:C
                matrix(:,c)=scalers{c}.inverse( matrix(:,c));
            end
        end
    end
end