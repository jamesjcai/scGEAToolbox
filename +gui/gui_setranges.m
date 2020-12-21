function gui_setranges(data)
% https://www.mathworks.com/matlabcentral/answers/143306-how-to-move-a-plotted-line-vertically-with-mouse-in-a-gui

fh=figure();
xdata=randn(300,1);
ydata=randn(300,1);
sh=scatter(xdata,ydata);
lh1=xline(0.3,'r-');
lh2=xline(0.7,'r-');
lh3=yline(0.3,'r-');
lh4=yline(0.7,'r-');

guidata(fh,[xdata ydata]);
%%
set(fh, ...
    'WindowButtonDownFcn',   @mouseDownCallback, ...
    'WindowButtonUpFcn',     @mouseUpCallback,   ...
    'WindowButtonMotionFcn', @mouseMotionCallback); 


function mouseDownCallback(figHandle,varargin)

    % get the handles structure
      xydata = guidata(figHandle);
%     lh1=handles{1}; lh2=handles{2};
%     lh3=handles{3}; lh4=handles{4};
    
    % get the position where the mouse button was pressed (not released)
    % within the GUI
    currentPoint = get(figHandle, 'CurrentPoint');
    x            = currentPoint(1,1);
    y            = currentPoint(1,2);
    
    % get the position of the axes within the GUI
    % allAxesInFigure = findall(figHandle,'type','axes');
    axes1 = get(figHandle, 'CurrentAxes');
    set(axes1,'Units','pixels');
    axesPos = get(axes1,'Position');
    minx    = axesPos(1);
    miny    = axesPos(2);
    maxx    = minx + axesPos(3);
    maxy    = miny + axesPos(4);
    
    % is the mouse down event within the axes?
    if x>=minx && x<=maxx && y>=miny && y<=maxy 

        % do we have graphics objects?
        % if isfield(handles,'plotHandles')
            
            % get the position of the mouse down event within the axes
            currentPoint = get(axes1, 'CurrentPoint');
            x            = currentPoint(2,1);
            y            = currentPoint(2,2);
            
            % b=findall(axes1.Children,'type','ConstantLine');
            % b(1).InterceptAxis
    if min(abs(x-axes1.XLim)) < min(abs(y-axes1.YLim)) 
            if abs(x-lh1.Value) < abs(x-lh2.Value)
                if ~isempty(lh1), delete(lh1); end
                lh1=xline(x,'r-');
            else
                if ~isempty(lh2), delete(lh2); end
                lh2=xline(x,'r-');
            end
    else
            if abs(y-lh3.Value) < abs(y-lh4.Value)
                if ~isempty(lh3), delete(lh3); end
                lh3=yline(y,'r-');
            else
                if ~isempty(lh4), delete(lh4); end
                lh4=yline(y,'r-');
            end
    end
    i=(xydata(:,1)>lh1.Value) & (xydata(:,1)<lh2.Value);
    j=(xydata(:,2)>lh3.Value) & (xydata(:,2)<lh4.Value);
    
    axes1.Title.String=sprintf('%d out of %d (%.2f%%)',...
        sum(i&j),length(i),100*sum(i&j)./length(i));
    end
    
end

function mouseUpCallback(figHandle,varargin)

    % get the handles structure
    handles = guidata(figHandle);
    
    if isfield(handles,'mouseIsDown')
        disp('mouseIsDown');
        if handles.mouseIsDown
            % reset all moving graphic fields
            handles.mouseIsDown     = false;
            handles.movingPlotHndle = [];
            handles.prevPoint       = [];
            
            % save the data
            guidata(figHandle,handles);
        end
    end
end

function mouseMotionCallback(figHandle,varargin)

    % get the handles structure
    handles = guidata(figHandle);
    
    if isfield(handles,'mouseIsDown')
        
        if handles.mouseIsDown
            currentPoint = get(handles.axes1, 'CurrentPoint');
            x            = currentPoint(2,1);
            y            = currentPoint(2,2);

            % compute the displacement from previous position to current
            xDiff = x - handles.prevPoint(1);
            yDiff = y - handles.prevPoint(2);

            % adjust this for the data corresponding to movingPlotHndle
            xData = get(handles.movingPlotHndle,'XData');
            yData = get(handles.movingPlotHndle,'YData');

            set(handles.movingPlotHndle,'YData',yData+yDiff,'XData',xData+xDiff);

            handles.prevPoint = [x y];
            
            % save the data
            guidata(figHandle,handles);
        end
    end
end

end