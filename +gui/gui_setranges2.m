function [idx,xr,yr]=gui_setranges2(x,y,xr,yr,txtx,txty)
% https://www.mathworks.com/matlabcentral/answers/143306-how-to-move-a-plotted-line-vertically-with-mouse-in-a-gui
if nargin<1, x=randn(300,1); end
if nargin<2, y=randn(300,1)*1000; end
if nargin<3, xr=[-.8 .8]; end
if nargin<4, yr=[-.8 .8]*1000; end
if nargin<5, txtx=''; end
if nargin<6, txty=''; end

fh=figure();
scatter(x,y);
lh1=xline(xr(1),'k:');
lh2=xline(xr(2),'r-');
lh3=yline(yr(1),'k:');
lh4=yline(yr(2),'r-');
idx=true(length(x),1);
guidata(fh,[x y]);
xlabel(txtx);
ylabel(txty);
set(fh,'WindowButtonDownFcn', @mouseDownCallback);

% ButtonH=uicontrol('Parent',fh,'Style','pushbutton',...
%     'String','Done','Units','normalized','Position',[0.0 0.5 0.4 0.2],'Visible','on');

       uicontrol(fh,'String',...
                  'Done',...
                  'Units','normalized',...
                  'Position',[0.45 0.0 0.2 0.056],...
                  'Callback', @i_CloseFig,...
                  'Tag','button');
              
       uicontrol(fh,'String',...
                  'Set Cutoffs...',...
                  'Units','normalized',...
                  'Position',[0.20 0.0 0.2 0.056],...
                  'Callback', @i_SetValues,...
                  'Tag','button');           

              
              waitfor(fh);


function i_CloseFig(figHandle,varargin)
    % delete(fh);
    closereq;
end

    function i_SetValues(figHandle,varargin)
            prompt = {'Library size cutoff:','Mt-ratio cutoff:'};
            answer = inputdlg(prompt,"",[1 35],...
                     {num2str(xr(2)),num2str(yr(2))});
            if isempty(answer), return; end
            try
                a1 = str2double(answer{1});
                a2 = str2double(answer{2});
                xr(2)=a1;
                yr(2)=a2;
                if ~isempty(lh2), delete(lh2); end
                lh2=xline(xr(2),'r-');
                if ~isempty(lh4), delete(lh4); end
                lh4=yline(yr(2),'r-');
                
                % xydata = guidata(fh);
                axes1 = get(fh,'CurrentAxes');
    i=(x>lh1.Value) & (x<lh2.Value);
    j=(y>lh3.Value) & (y<lh4.Value);
    idx=i&j;
    xr=[lh1.Value lh2.Value];
    yr=[lh3.Value lh4.Value];
    axes1.Title.String=sprintf('Inclusion: %d out of %d (%.2f%%)\nExclusion: %d out of %d (%.2f%%)',...
        sum(i&j),length(i),100*sum(idx)./length(i),...
        length(i)-sum(i&j),length(i),100*(length(i)-sum(idx))./length(i));
            catch
                errordlg('Wrong inputs')
                return;
            end            
    end

function mouseDownCallback(figHandle,varargin)

    % get the handles structure
      xydata = guidata(figHandle);
%     lh1=handles{1}; lh2=handles{2};
%     lh3=handles{3}; lh4=handles{4};
    
    % get the position where the mouse button was pressed (not released)
    % within the GUI
    currentPoint = get(figHandle, 'CurrentPoint');
    x1            = currentPoint(1,1);
    y1            = currentPoint(1,2);
    
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
    if x1>=minx && x1<=maxx && y1>=miny && y1<=maxy 

        % do we have graphics objects?
        % if isfield(handles,'plotHandles')
            
            % get the position of the mouse down event within the axes
            currentPoint = get(axes1, 'CurrentPoint');
            xx            = currentPoint(2,1);
            yy            = currentPoint(2,2);
            
            % b=findall(axes1.Children,'type','ConstantLine');
            % b(1).InterceptAxis
            
    if min(abs((xx-axes1.XLim)./diff(axes1.XLim))) < ...
       min(abs((yy-axes1.YLim)./diff(axes1.YLim))) 
            if abs(xx-lh1.Value) < abs(xx-lh2.Value)
                if ~isempty(lh1), delete(lh1); end
                %lh1=xline(xx,'r-');
                lh1=xline(0,'k:');
            else
                if ~isempty(lh2), delete(lh2); end
                lh2=xline(xx,'r-');
            end
    else
            if abs(yy-lh3.Value) < abs(yy-lh4.Value)
                if ~isempty(lh3), delete(lh3); end
                %lh3=yline(yy,'r-');
                lh3=yline(0,'k:');
            else
                if ~isempty(lh4), delete(lh4); end
                lh4=yline(yy,'r-');
            end
    end
    i=(xydata(:,1)>lh1.Value) & (xydata(:,1)<lh2.Value);
    j=(xydata(:,2)>lh3.Value) & (xydata(:,2)<lh4.Value);
    idx=i&j;
    xr=[lh1.Value lh2.Value];
    yr=[lh3.Value lh4.Value];
    axes1.Title.String=sprintf('Inclusion: %d out of %d (%.2f%%)\nExclusion: %d out of %d (%.2f%%)',...
        sum(i&j),length(i),100*sum(idx)./length(i),...
        length(i)-sum(i&j),length(i),100*(length(i)-sum(idx))./length(i));
    end
end

end