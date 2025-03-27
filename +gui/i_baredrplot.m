function i_baredrplot(ax, c, t, parentfig)

if nargin < 4, parentfig = []; end
if nargin < 3, t = 'tSNE'; end
if nargin < 2, c = []; end

% if size(s, 2) >= 3
%     x = s(:, 1);
%     y = s(:, 2);
%     z = s(:, 3);
%     is2d = false;
% else
%     x = s(:, 1);
%     y = s(:, 2);
%     z = zeros(size(x));
%     is2d = true;
% end
isAxesHandle = isa(ax, 'matlab.graphics.axis.Axes'); %isgraphics(s, 'axes');
if ~isAxesHandle && isempty(c), error('Empty handle.'); end

hx = gui.myFigure(parentfig);
% hx = gui.myFigure;
hFig = hx.FigHandle;

%if ~isempty(parentfig)
%    hFig.Position = parentfig.Position;
%end

if isAxesHandle
    if gui.i_isuifig(hFig)
        copyobj(ax.Children, hx.AxHandle);
        hAx = hx.AxHandle;
    else
        hAx = copyobj(ax, hFig);
    end
else
    if gui.i_isuifig(hFig)
        hAx = hx.AxHandle;
    else
        hAx = axes('Parent', hFig, 'Visible', 'off');
    end
    %h1 = gui.i_gscatter3(ax, c, 1, 1, hAx);  
end


set(hAx, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');
grid(hAx, 'off');

% Get axis limits and figure position
% ax = gca;
% ax = hAx;
xLimits = hAx.XLim;
yLimits = hAx.YLim;
zLimits = hAx.ZLim;

%bx=gca;
%assert(isequal(bx, hAx));

is3d1 = isprop(hAx, 'ZLim');
h = get(hAx, 'Children');
is3d2 = false;
if isscalar(h) 
    if isprop(h, 'ZData')  % ismember('ZData', properties(h))
        is3d2 = any(arrayfun(@(x) ~isempty(get(x, 'ZData')), h));
    end
else    
    for k=1:numel(h)
        if isa(h(k), 'matlab.graphics.chart.primitive.Scatter')
            is3d2 = any(arrayfun(@(x) ~isempty(get(x, 'ZData')), h(k)));
        end
    end
end
is3d = is3d1 & is3d2;

hold(hAx, 'on');

if is3d    % ======================================== 3D
    % Turn off the default axis display
    %set(gca, 'Visible', 'off');
 
    a = xLimits(1); b = yLimits(1); c = zLimits(1);
    la = xLimits(2)-a;
    lb = yLimits(2)-b;
    lc = zLimits(2)-c;

    % Draw custom arrows as axes using quiver3
    quiver3(hAx, a, b, c, la/5, 0, 0, 'k', 'LineWidth', 1); % X-axis
    quiver3(hAx, a, b, c, 0, lb/5, 0, 'k', 'LineWidth', 1); % Y-axis
    quiver3(hAx, a, b, c, 0, 0, lc/5, 'k', 'LineWidth', 1); % Z-axis

%    axis_length = 20;
%    quiver3(0, 0, 0, axis_length, 0, 0, 'k', 'LineWidth', 1); % X-axis
%    quiver3(0, 0, 0, 0, axis_length, 0, 'k', 'LineWidth', 1); % Y-axis
%    quiver3(0, 0, 0, 0, 0, axis_length, 'k', 'LineWidth', 1); % Z-axis
    txt1 = sprintf('%s\\_1', t);
    txt2 = sprintf('%s\\_2', t);
    txt3 = sprintf('%s\\_3', t);
    
    % Label each arrow for clarity
    text(hAx, a+la/5, b, c, txt1);
    text(hAx, a, b+lb/5, c, txt2);
    text(hAx, a, b, c+lc/5, txt3);
    view(hAx, 3);
else          % ======================================== 2D
    %disp('2D')
    hAx.Units = "pixels";
    r = hAx.Position(3)/hAx.Position(4);
    
    hAx.Units = "normalized";
    axPos = hAx.Position;
    
    % Convert data limits to figure normalized units
    % xArrowPos = [0, 1]; % normalized from left to right of the axes
    % yArrowPos = [0, 1]; % normalized from bottom to top of the axes
    
    % Draw x-axis arrow
    annotation(hFig, 'arrow', ...
        [axPos(1), axPos(1) + axPos(3)/7], ... % x positions
        [axPos(2), axPos(2)], ...            % y positions
        'Color', 'k', 'LineWidth', .5);
    
    % Draw y-axis arrow
    annotation(hFig, 'arrow', ...
        [axPos(1), axPos(1)], ...            % x positions
        [axPos(2), axPos(2) + r*(axPos(4)/7)], ... % y positions
        'Color', 'k', 'LineWidth', .5);
        
    txt1 = sprintf('%s\\_1', t);
    txt2 = sprintf('%s\\_2', t);

    textOpts = struct();
    textOpts.HorizontalAlignment = 'center';
    textOpts.VerticalAlignment = 'middle';
    textOpts.FontSize = 10;
    textOpts.FontWeight = 'normal';

    [~, b] = measureText(txt1, textOpts, hAx);
    text(hAx, xLimits(1), yLimits(1) - 2*b, txt1);
    
    [~, b] = measureText(txt2, textOpts, hAx);
    text(hAx, xLimits(1)-3*b, yLimits(1), txt2,'Rotation',90);
    view(hAx, 2);
end

hold(hAx,'off');
title(hAx, '')
subtitle(hAx, '')
xlim(hAx, xLimits);
ylim(hAx, yLimits);
hx.show(parentfig);


 function [width, height] = measureText(txt, textOpts, ax)
    % if(nargin < 3)
    %    ax = gca();
    % end
    % if nargin < 2
    %     textOpts = struct();
    %     textOpts.HorizontalAlignment = 'center';
    %     textOpts.VerticalAlignment = 'middle';
    %     textOpts.FontSize = 10;
    %     textOpts.FontWeight = 'normal';
    % end
    hTest = text(ax, 0, 0, txt, textOpts);
    textExt = get(hTest, 'Extent');
    delete(hTest);
    height = textExt(4)/3;    %Height
    width = textExt(3)/3;     %Width
 end

end
