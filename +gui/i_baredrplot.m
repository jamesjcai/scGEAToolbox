function [h1] = i_baredrplot(s, c, t, parentfig)

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
isAxesHandle = isa(s, 'matlab.graphics.axis.Axes'); %isgraphics(s, 'axes');
if ~isAxesHandle && isempty(c), error('xxx'); end


hFig = figure("Visible","off");

if isAxesHandle
    hAx = copyobj(s, hFig);
else
    hAx = axes('Parent', hFig, 'Visible', 'off');
    h1 = gui.i_gscatter3(s, c, 1, 1, hAx);    
end
grid off
hold on;

set(hAx, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');
% Get axis limits and figure position
ax = gca;
xLimits = ax.XLim;
yLimits = ax.YLim;
zLimits = ax.ZLim;

assert(isequal(ax, hAx));

is3d1 = isprop(hAx, 'ZLim');
h = get(hAx, 'Children');
is3d2 = any(arrayfun(@(x) ~isempty(get(x, 'ZData')), h));
is3d = is3d1 & is3d2;

if is3d    % ======================================== 3D
    hold on;
    % Turn off the default axis display
    %set(gca, 'Visible', 'off');
 
    a = xLimits(1); b = yLimits(1); c = zLimits(1);
    la = xLimits(2)-a;
    lb = yLimits(2)-b;
    lc = zLimits(2)-c;

    % Draw custom arrows as axes using quiver3
    quiver3(a, b, c, la/5, 0, 0, 'k', 'LineWidth', 1); % X-axis
    quiver3(a, b, c, 0, lb/5, 0, 'k', 'LineWidth', 1); % Y-axis
    quiver3(a, b, c, 0, 0, lc/5, 'k', 'LineWidth', 1); % Z-axis

    axis_length = 20;
%    quiver3(0, 0, 0, axis_length, 0, 0, 'k', 'LineWidth', 1); % X-axis
%    quiver3(0, 0, 0, 0, axis_length, 0, 'k', 'LineWidth', 1); % Y-axis
%    quiver3(0, 0, 0, 0, 0, axis_length, 'k', 'LineWidth', 1); % Z-axis
    txt1 = sprintf('%s\\_1', t);
    txt2 = sprintf('%s\\_2', t);
    txt3 = sprintf('%s\\_3', t);
    
    % Label each arrow for clarity
    text(a+la/5, b, c, txt1);
    text(a, b+lb/5, c, txt2);
    text(a, b, c+lc/5, txt3);
    hold off;
    

else          % ======================================== 2D
    ax.Units = "pixels";
    r = ax.Position(3)/ax.Position(4);
    
    ax.Units = "normalized";
    axPos = ax.Position;
    
    % Convert data limits to figure normalized units
    % xArrowPos = [0, 1]; % normalized from left to right of the axes
    % yArrowPos = [0, 1]; % normalized from bottom to top of the axes
    
    % Draw x-axis arrow
    annotation('arrow', ...
        [axPos(1), axPos(1) + axPos(3)/7], ... % x positions
        [axPos(2), axPos(2)], ...            % y positions
        'Color', 'k', 'LineWidth', .5);
    
    % Draw y-axis arrow
    annotation('arrow', ...
        [axPos(1), axPos(1)], ...            % x positions
        [axPos(2), axPos(2) + r*(axPos(4)/7)], ... % y positions
        'Color', 'k', 'LineWidth', .5);
    %axis off
    
    txt1 = sprintf('%s\\_1', t);
    txt2 = sprintf('%s\\_2', t);
    
    [~, b] = measureText(txt1);
    text(xLimits(1), yLimits(1) - 2*b, txt1);
    
    [~, b] = measureText(txt2);
    text(xLimits(1) - 2*b, yLimits(1), txt2,'Rotation',90);
    title('')
    subtitle('')
    hold off;
end



gui.i_movegui2parent(hFig, parentfig);
hFig.Visible = "on";



 function [width, height] = measureText(txt, textOpts, axis)
    if(nargin < 3)
       axis = gca(); 
    end
    if nargin < 2
        textOpts = struct();
        textOpts.HorizontalAlignment = 'center';
        textOpts.VerticalAlignment = 'middle';
        textOpts.FontSize = 10;
        textOpts.FontWeight = 'normal';
    end
    hTest = text(axis, 0, 0, txt, textOpts);
    textExt = get(hTest, 'Extent');
    delete(hTest);
    height = textExt(4)/3;    %Height
    width = textExt(3)/3;     %Width
end

end