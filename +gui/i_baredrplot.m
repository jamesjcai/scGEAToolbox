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

% Turn off default x and y axes
set(gca, 'XColor', 'none', 'YColor', 'none');
set(gca, 'ZColor', 'none');
% Get axis limits and figure position
ax = gca;

xLimits = ax.XLim;
yLimits = ax.YLim;
ax.Units = "pixels";
r = ax.Position(3)/ax.Position(4);

ax.Units = "normalized";
axPos = ax.Position;

% Convert data limits to figure normalized units
xArrowPos = [0, 1]; % normalized from left to right of the axes
yArrowPos = [0, 1]; % normalized from bottom to top of the axes

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