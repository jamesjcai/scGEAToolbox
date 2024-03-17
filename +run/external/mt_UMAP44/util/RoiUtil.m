classdef RoiUtil < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer: Connor Meehan <connor.gw.meehan@gmail.com>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(Constant)
        V='9.6';
        r2020a='9.8';
        r2021b='9.11';
        POLYGON='impoly';
        RECTANGLE='imrect';
        ELLIPSE='imellipse';
        ELLIPSE_NEW='images.roi.Ellipse';
        RECTANGLE_NEW='images.roi.Rectangle';
        POLYGON_NEW='images.roi.Polygon';
        NEW_COLOR=[.41 .86 .9];
        EDIT_COLOR=[.91 .91 .72];
    end
    
    methods(Static)
        function ok=CanDoNew
            if ~verLessThan('matlab', RoiUtil.r2021b)
                ok=true;
                return;
            end
            app=BasicMap.Global;
            if isempty(app.newRoiBroken)
                RoiUtil.InitRoi(app); % will fail once
            end
            ok=~verLessThan('matlab',  RoiUtil.V) && ~app.newRoiBroken;
        end
        
        function [roi, str]=NewForXy(ax, xy, tolerance, cbMoved)
            if nargin<4
                cbMoved=[];
                if nargin<3
                    tolerance=1;
                end
            end
            xy=xy(boundary(xy(:,1), xy(:,2), .1241), :);
            xy=double(edu.stanford.facs.swing.CurveSimplifier.ToArray(...
                xy, tolerance));
            [roi, str]=RoiUtil.NewPolygon(ax, xy, cbMoved);
        end
        
        function oldClr=SetColor(roi, clr)
            oldRoi=~RoiUtil.CanDoNew;
            if ~oldRoi
                oldClr=roi.Color;
                roi.Color=clr;
            else
                oldClr=roi.getColor;
                roi.setColor(clr);
            end
        end
        
        function oldLbl=SetLabel(roi, lbl)
            oldRoi=~RoiUtil.CanDoNew;
            try
                if ~oldRoi
                    oldLbl=roi.Label;
                    roi.Label=lbl;
                end
            catch
            end
        end
        
        
        function [roi, str]=NewPolygon(ax, xy, cbMoved)
            str=MatBasics.XyToString(xy);            
            oldRoi=~RoiUtil.CanDoNew;
            if ~oldRoi
                roi=drawpolygon(ax, ...
                    'Position', xy,...
                    'ContextMenu', [], 'SelectedColor', ...
                    RoiUtil.NEW_COLOR);
            else
                roi=impoly(ax, xy);
                roi.setColor(RoiUtil.NEW_COLOR);
            end
            if ~isempty(cbMoved)
                if oldRoi
                    roi.addNewPositionCallback(...
                        @(pos)RoiUtil.OldRoiMoved(cbMoved, roi));
                else
                    addlistener(roi, 'ROIMoved', ...
                        @(src,evt)RoiUtil.Moved(cbMoved, src));
                    
                end
                feval(cbMoved,roi);
            end
        end
        
        function roi=New(ax, type, cbMoved, cbMoving)
            oldRoi=~RoiUtil.CanDoNew;
            if strcmp(type, RoiUtil.RECTANGLE)
                if ~oldRoi
                    roi=drawrectangle(ax, 'ContextMenu', [],...
                        'Rotatable', true, 'SelectedColor', ...
                        RoiUtil.NEW_COLOR);
                else
                    roi=imrect(ax);
                    RoiUtil.SetColor(roi, RoiUtil.NEW_COLOR);
                end
            elseif strcmp(type, RoiUtil.ELLIPSE)
                if ~oldRoi %&& ~RoiUtil.MustUseOldRoi(ax)
                    roi=drawellipse(ax, 'RotationAngle', 0, ...
                        'ContextMenu', [], 'SelectedColor', ...
                        RoiUtil.NEW_COLOR);
                else
                    roi=imellipse(ax);
                    RoiUtil.SetColor(roi, RoiUtil.NEW_COLOR);
                end
            else
                if ~oldRoi
                    roi=drawpolygon(ax, ...
                        'ContextMenu', [], 'SelectedColor', ...
                        RoiUtil.NEW_COLOR);
                else
                    roi=impoly(ax);
                    RoiUtil.SetColor(roi, RoiUtil.NEW_COLOR);
                end
            end
            if nargin>2                
                try
                    if oldRoi
                        roi.addNewPositionCallback(...
                            @(pos)RoiUtil.OldRoiMoved(cbMoved, roi));
                    else
                        if nargin>3
                            addlistener(roi, 'MovingROI', ...
                                @(src,evt)RoiUtil.Moving(cbMoving, src));
                        end
                        addlistener(roi, 'ROIMoved', ...
                            @(src,evt)RoiUtil.Moved(cbMoved, src));
                    end
                    feval(cbMoved,roi);
                catch ex
                    ex.getReport
                    roi=[];
                end
            end
        end
        
        function must=MustUseOldRoi(ax)
            yl=ylim(ax);
            xl=xlim(ax);
            yRange=abs(yl(2) - yl(1));
            xRange=abs(xl(2)-xl(1));
            dif=yRange/xRange;
            must=dif<.66 || dif>1.33;
            if ~must
                 % weird performance issue with ellipse ROI
                must=xRange>5 || yRange>5; 
            end
        end
        
        
        function OldRoiMoved(cb, roi)
            feval(cb, roi);
        end
        
        function Moving(cb, roi)
            feval(cb, roi);
        end

        function Moved(cb, roi)
            feval(cb, roi);
        end

        function Clicked(cb, roi, evt)
            feval(cb, roi, evt);
        end

        function ok=IsNewRoi(roi)
            if ischar(roi)
                ok=strcmp(roi, 'images.roi.Ellipse') ||...
                    strcmp(roi, 'images.roi.Rectangle') ||...
                    strcmp(roi, 'images.roi.Polygon');
            else
                ok=isa(roi, 'images.roi.Ellipse') ||...
                    isa(roi, 'images.roi.Rectangle') ||...
                    isa(roi, 'images.roi.Polygon');
            end
        end

        function yes=IsHandle(roi)
            if RoiUtil.IsNewRoi(roi)
                yes=ishandle(roi);
            else
                try
                    roi.getPosition;
                    yes=true;
                catch
                    yes=false;
                end
            end
        end
        
        function position=Position(roi)
            if isempty(roi)
                position=[];
                return;
            end
            if ~RoiUtil.IsNewRoi(roi)
                position=roi.getPosition();
            else
                if isa(roi, 'images.roi.Ellipse')
                    c=get(roi, 'Center');
                    sa=get(roi, 'SemiAxes');
                    position=[c(1)-sa(1) c(2)-sa(2) sa(1)*2 sa(2)*2];
                    if roi.RotationAngle ~=0
                        position(end+1)=roi.RotationAngle;
                    end
                else
                    position=get(roi, 'Position');
                    if isa(roi, 'images.roi.Rectangle')
                        if roi.RotationAngle ~=0
                            position(end+1)=roi.RotationAngle;
                        end
                    end
                end
            end
        end
        
                
        function rows=GetRows(roi, data2D)
            newRoi=RoiUtil.IsNewRoi(roi);
            if newRoi 
                rows=roi.inROI(data2D(:,1), data2D(:,2));
                return;
            end
            pos=RoiUtil.Position(roi);
            typeT=class(roi);
            if strcmpi(typeT,RoiUtil.ELLIPSE)
                rows=RoiUtil.InEllipseUnrotated(data2D(:,1), data2D(:,2), pos);
            elseif strcmpi(typeT, RoiUtil.RECTANGLE)
                rows=RoiUtil.InRectUnrotated(data2D(:,1), data2D(:,2), pos);
            elseif strcmp(typeT, RoiUtil.POLYGON)
                rows=inpolygon(data2D(:,1), data2D(:,2),pos(:,1),pos(:,2));
            else
                rows=[];
            end
        end
        

        function inside=InEllipseUnrotated(X, Y, pos)
            xmin = pos(1);
            ymin = pos(2);
            width = pos(3);
            height = pos(4);
            a = width/2;
            b = height/2;
            center = [xmin+a, ymin + b];
            inside = (X - center(1)*ones(size(X))).^2./a^2 + ...
                (Y - center(2)*ones(size(Y))).^2./b^2 <= ones(size(X));
        end
        
        function inside=InRectUnrotated(X, Y, pos)
            xmin = pos(1);
            ymin = pos(2);
            width = pos(3);
            height = pos(4);
            a = width/2;
            b = height/2;
            center = [xmin+a, ymin + b];
            inside = abs(X - center(1)*ones(size(X))) <= a*ones(size(X))...
                & abs(Y - center(2)*ones(size(Y))) <= b*ones(size(Y));
        end
        
        function ok=IsHandleOk(roi)
            try
                ok=ishandle(roi.Parent);
                if isempty(ok) %odd .... does NOT look ok
                    ok=false;
                end
            catch
                ok=false;
            end
        end
        
        function InitRoi(app)
            app.oldRoi=verLessThan('matlab',  RoiUtil.V);
            app.newRoiBroken=[];
            if ~app.oldRoi
                MatBasics.RunLater(@(h,e)testROIStep1, .31);
            end
            
            function testROIStep1
                R=RoiUtil.GetRectRoi(app);
                fig=get(get(R, 'parent'), 'parent'); 
                MatBasics.RunLater(@(h,e)testROIStep2(fig), .31);
                R.inROI([1 2], [2 1]);
                app.newRoiBroken=false;                
            end
            
            function testROIStep2(fig)
                if isempty(app.newRoiBroken)
                    app.newRoiBroken=true;
                    app.oldRoi=true;
                    warning(['Reverting to pre-2018b regions of interest '...
                        'since inROI() in images.roi.* does not work!']);
                end
                delete(fig);
            end
            
        end
        
        function [R, ax, fig]=NewRect(pos, ax, newFigIfPossible)
            if nargin<3
                newFigIfPossible=true;
                if nargin<2
                    ax=[];
                end
            end
            app=BasicMap.Global;
            if isempty(ax)
                fig=figure('visible', 'off', ...
                    'name', 'dummyRect');
                ax=Gui.Axes(fig);
            else
                fig=[];
            end
            if newFigIfPossible && ~app.oldRoi                
                R=drawrectangle(ax, ...
                    'position', pos, ...
                    'Rotatable', true, ...
                    'ContextMenu', [],...
                    'InteractionsAllowed', 'none');
            else
                R=imrect(ax, pos);
            end
        end
        
        function R=GetRectRoi(app)
            if ~app.oldRoi
                if isempty(app.rectRoi)
                    app.rectRoi=drawrectangle(...
                        Gui.Axes(figure('visible', 'off', ...
                        'name', 'dummyRect')), ...
                        'position', [10 20 30 90], ...
                        'Rotatable', true, ...
                        'ContextMenu', [],...
                        'InteractionsAllowed', 'none');
                end
            end
            R=app.rectRoi;
        end
        
        function R=GetEllipseRoi(app)
            if isempty(app.ellipseRoi)
                if ~app.oldRoi
                    app.ellipseRoi=drawellipse(...
                        Gui.Axes(figure('visible', 'off', ...
                        'name', 'dummyEllipse')), ...
                        'Center', [10 10], ...
                        'SemiAxes', [5 5], ...
                        'InteractionsAllowed', 'none');
                end
            end
            R=app.ellipseRoi;
        end

        function position=ToPositionFromFlowJoEllipse(edge, scalers)
            [center, semiAxes, rotationAngle]=RoiUtil.ToRoiEllipse(edge);
            binPosition=[center(1)-semiAxes(1) center(1)+semiAxes(1);...
                center(2)-semiAxes(2) center(2)+semiAxes(2)];
            position(1,:)=scalers{1}.inverseFromBins(binPosition(1, :), 256);
            position(2,:)=scalers{2}.inverseFromBins(binPosition(2, :), 256);
            position=[position(1,1) position(2,1) ...
                abs(position(1,2)-position(1,1)) ...
                abs(position(2,2)-position(2,1))...
                rotationAngle];
        end

        function [center, semiAxes, rotationAngle]=...
                ToRoiEllipse(edge)
            %r2018b or later images.roi.Ellipse
            originVertex = edge(1,:); oppositeOriginVertex = edge(2,:);
            otherVertex = edge(3,:); otherVertex2 = edge(4,:);
            
            center = (originVertex+oppositeOriginVertex)/2;

            semiAxes = [norm(originVertex - oppositeOriginVertex)/2 ...
                        norm(otherVertex - otherVertex2)/2];

            rotationAngle = MatBasics.getPositiveClockwiseAngleD(originVertex,oppositeOriginVertex);
        end

        function [foci, edge]=...
                ToFlowJoEllipseFromPosition(pos, scalers)
            [xMin, xMax, yMin, yMax]=RoiUtil.ToFlowJoRect(pos);
            binPos=[xMin xMax;yMin yMax];
            if nargin<2
                binPos=SuhScaler.ScaledToBins(binPos, 256);
            else
                binPos(:,1)=scalers{1}.unscaledToBins(binPos(:,1), 256);
                binPos(:,2)=scalers{2}.unscaledToBins(binPos(:,2), 256);
            end
            W=binPos(1,2)-binPos(1,1);
            H=binPos(2,2)-binPos(2,1);
            semiAxes=[W/2 H/2];
            center=[binPos(1,1)+semiAxes(1) binPos(2,1)+semiAxes(2)];
            if length(pos)>4
                rotationAngle=pos(5);
            else
                rotationAngle=0;
            end
            [foci, edge]=RoiUtil.ToFlowJoEllipse(...
                center, semiAxes, rotationAngle);
            edge=round(edge);
        end

        function [foci, edge]=...
                ToFlowJoEllipse(center, semiAxes, rotationAngle)
            originSemiAxis = semiAxes(1); otherSemiAxis = semiAxes(2);
            majorSemiAxis = max(semiAxes); minorSemiAxis = min(semiAxes);
            counterClockwiseAngle = 360-rotationAngle;

            vertexAngles = counterClockwiseAngle+(0:90:270);
            vertexDirections = [cosd(vertexAngles)' sind(vertexAngles)'];

            flowJoVertexOrder = [3; 1; 4; 2]; 
         
            originVertex = center + originSemiAxis*vertexDirections(flowJoVertexOrder(1),:);
            oppositeOriginVertex = center + originSemiAxis*vertexDirections(flowJoVertexOrder(2),:);
            otherVertices = center + otherSemiAxis*vertexDirections(flowJoVertexOrder(3:4),:);

            [~,I1] = sort(otherVertices(:,2));      %For the third and fourth vertices, FlowJo appears to 
            otherVertices = otherVertices(I1,:);    %list the lower vertex first.

            edge = [originVertex; oppositeOriginVertex; otherVertices];

            focusOffset = sqrt(majorSemiAxis^2 - minorSemiAxis^2);
            fociAlongOriginSemiAxis = (originSemiAxis == majorSemiAxis);

            if fociAlongOriginSemiAxis
                foci = center + focusOffset*vertexDirections(flowJoVertexOrder(1:2),:);
            else
                foci = center + focusOffset*vertexDirections(flowJoVertexOrder(3:4),:);
            end

            [~,I2] = sort(foci(:,1));  %FlowJo appears to list the leftmost
            foci = foci(I2,:);         %focus first.
        end

        function [xMin, xMax, yMin, yMax]=ToFlowJoRect(roiPosition)
            xMin=roiPosition(1);
            yMin=roiPosition(2);
            xMax=roiPosition(1)+roiPosition(3);
            yMax=roiPosition(2)+roiPosition(4);
        end
        
        function yes=IsEllipse(roi)
            type=class(roi);
            yes=strcmp(type, RoiUtil.ELLIPSE) || ...
                strcmp(type, RoiUtil.ELLIPSE_NEW);
        end
    end
    
    properties
        position;
        label;
        roi;
        type;
        newRoi;%r2018b or later
        cbMoved;
        isHistogram;
        is1DRectangle;
    end
    
    methods
        function select(this, yes)
            if nargin<2
                yes=true;
            end
            if this.newRoi
                this.roi.Selected=yes;
                if yes
                    bringToFront(this.roi)
                end
            end
        end
        
        function ok=isSelected(this)
            try
                if this.newRoi
                    ok=this.roi.Selected;
                else
                    ok=false;
                end
            catch
                ok=false;
            end
        end
        

        function setLabel(this, label)
            this.label=label;
            if this.newRoi
                if this.IsHandleOk(this.roi)
                    this.roi.Label=label;
                end
            end
        end
        
        function this=RoiUtil(type, pos, isHistogram)
            if nargin<3
                isHistogram=false;
            end
            if isa(type, 'RoiUtil')
                prior=type;
                type=prior.type;
                if nargin<2
                    pos=prior.position;
                end
            else
                prior=[];
            end
            this.newRoi=RoiUtil.CanDoNew;
            this.type=type;
            this.isHistogram=isHistogram;
            if length(pos)==2
                if strcmp(type, RoiUtil.RECTANGLE)
                    pos=[pos(1) pos(1) pos(2) pos(2)];
                    this.is1DRectangle=~ isHistogram;
                end
            end
            this.position=pos;
            if this.newRoi
                if strcmp(type, RoiUtil.RECTANGLE)
                    this.roi=images.roi.Rectangle('Position', pos);
                elseif strcmp(type, RoiUtil.ELLIPSE)
                    semiAxes=[pos(3)/2 pos(4)/2];
                    center=[pos(1)+semiAxes(1) pos(2)+semiAxes(2)];
                    if length(pos)==5
                        this.roi=images.roi.Ellipse('RotationAngle', pos(5),...
                            'Center', center, 'SemiAxes', semiAxes);
                    else
                        this.roi=images.roi.Ellipse(...
                            'Center', center, 'SemiAxes', semiAxes);
                    end
                else
                    this.roi=images.roi.Polygon('Position', pos);
                end                
                if  ~isempty(prior)
                    this.roi.Label=prior.roi.Label;
                    this.label=prior.label;
                    this.roi.UserData=prior.roi.UserData;
                    this.roi.Color=prior.roi.Color;
                    this.roi.SelectedColor=prior.roi.SelectedColor;
                    this.setMovedCallback(prior.cbMoved);
                end
                try
                    this.roi.LabelAlpha=.2;
                    this.roi.LabelTextColor=[.45 .1 .68];
                    this.roi.MarkerSize=6;
                    this.roi.SelectedColor=RoiUtil.EDIT_COLOR;
                catch
                end
            end
        end
        
        function H=setParent(this, parent)
            if ~this.newRoi
                if strcmp(this.type, RoiUtil.RECTANGLE)
                    H=imrect(parent, this.position);
                elseif strcmp(this.type, RoiUtil.ELLIPSE)
                    H=imellipse(parent, this.position(1:4));
                else
                    H=impoly(parent, this.position);
                end
                this.roi=H;
            else
                try
                    H=this.roi;
                    H.Parent=parent;
                catch
                end
            end
        end
        
        function ok=isHandleOk(this)
            if ~this.newRoi
                ok=true;
            else
                try
                    ok=ishandle(this.roi.Parent);
                    if isempty(ok) %odd .... does NOT look ok
                        ok=false;
                    end
                catch
                    ok=false;
                end
            end
        end
        
        function pos=getPosition(this)
            pos=RoiUtil.Position(this.roi);
        end
        
        function setPosition(this, pos, setOnlyRoiObject)
            if nargin<3
                setOnlyRoiObject=false;
            end
            if ~setOnlyRoiObject
                this.position=pos;
            end
            if this.newRoi
                if strcmp(this.type, RoiUtil.RECTANGLE)
                    this.roi.Position=pos(1:4);
                    if length(pos)==5
                        this.roi.RotationAngle=pos(5);
                    end
                elseif strcmp(this.type, RoiUtil.ELLIPSE)
                    center=[pos(1)+pos(3)/2 pos(2)+pos(4)/2];
                    semiAxes=[pos(3)/2 pos(4)/2];
                    if length(pos)==5
                        this.roi.RotationAngle=pos(5);
                    end
                    this.roi.Center=center;
                    this.roi.SemiAxes=semiAxes;
                else
                    this.roi.Position=pos;
                end
            end
        end
        
        function rows=inROI(this, XY)
            %handle use of rectangle on 1D
            if isempty(XY{2})
                XY{2}=XY{1};
            elseif isempty(XY{1})
                XY{1}=XY{2};
            end
            if this.newRoi
                if strcmp(this.type, RoiUtil.POLYGON)
                    rows=inpolygon(XY{1}, XY{2},...
                        this.position(:,1), this.position(:,2));
                else
                    rows=this.roi.inROI(XY{1}, XY{2});
                end
                return;
            end            
            if strcmpi(this.type,RoiUtil.ELLIPSE)
                rows=RoiUtil.InEllipseUnrotated(XY{1}, XY{2}, ...
                    this.position);
            elseif strcmpi(this.type, RoiUtil.RECTANGLE)
                rows=RoiUtil.InRectUnrotated(XY{1}, XY{2}, ...
                    this.position);
            elseif strcmp(this.type, RoiUtil.POLYGON)
                rows=inpolygon(XY{1}, XY{2},...
                    this.position(:,1), this.position(:,2));
            else
                rows=[];
            end
        end
        
        function setMovedCallback(this, cbMoved)
            if ~isempty(cbMoved)
                this.cbMoved=cbMoved;
                if ~this.newRoi
                    this.roi.addNewPositionCallback(...
                        @(pos)RoiUtil.OldRoiMoved(cbMoved, ...
                        this.roi));
                else
                    addlistener(this.roi, 'ROIMoved', ...
                        @(src,evt)RoiUtil.Moved(cbMoved, this.roi));
                end
            end
        end
        
        function setClickedCallback(this, cbClicked)
            if ~isempty(cbClicked)
                this.cbMoved=cbClicked;
                if ~this.newRoi
                    disp('No clicked callback on ROI BEFORE r2018b');
                else
                    addlistener(this.roi, 'ROIClicked', ...
                        @(src,evt)RoiUtil.Clicked(cbClicked, src, evt));
                end
            end
        end
        
    end
end