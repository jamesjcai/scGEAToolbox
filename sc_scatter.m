function sc_scatter(varargin)
%SC_SCATTER Single-cell Gene Expression Analysis Toolbox Graphical User Interface.
%
% sc_scatter has been removed. Use scgeatool instead.
%
%   See also SCGEATOOL.

% Make an "error-dialog-like" dialog that can launch scgeatool if desired.

dlgWidth = 300;
dlgHeight = 150;

msgWidth = 190;
msgHeight = 50;
msgPosX = 60;
msgPosY = 75;
IconPosX = 12;
IconPosY = 62;
IconWidth = 24;
IconHeight = 24;

btnBaseWidth = 120;
btnBaseHeight = 25;
btnBasePosX = round(dlgWidth/2-btnBaseWidth/2);
btnBasePosY = 30;

d = dialog('Position', [400, 400, dlgWidth, dlgHeight], 'Name', '', ...
    'WindowStyle', 'modal', 'Tag', 'DeprecationDialog', 'Visible', false);

btn1 = uicontrol('Parent', d, ...
    'Position', [btnBasePosX, btnBasePosY, btnBaseWidth, btnBaseHeight], ...
    'String', 'Run SCGEATOOL', ...
    'Callback', @(~, ~)openOptimizeTaskCallback(d), ...
    'Tag', 'OpenOptimize');
buttonFitToText(btn1, 10);


% Resize the dialog box if needed
buffer = 10;
btnEdge = btn1.Position(1) + btn1.Position(3);
if btnEdge > d.Position(3) - buffer
    d.Position(3) = btnEdge + buffer;
end

uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [msgPosX, msgPosY, msgWidth, msgHeight], ...
    'String', 'SC_SCATTER has been renamed to SCGEATOOL.');


% Make the warning icon
% NOTE: much of what's here is cribbed from msgbox source to get the
% warning icon.

iconPos = [IconPosX, IconPosY, IconWidth, IconHeight];
IconAxes = axes( ...
    'Parent', d, ...
    'Units', 'points', ...
    'Position', iconPos, ...
    'Tag', 'IconAxes' ...
    );

[iconData, alphaData] = matlab.ui.internal.dialog.DialogUtils.imreadDefaultIcon('error');
Img = image('CData', iconData, 'Parent', IconAxes);
if ~isempty(alphaData)
    set(Img, 'AlphaData', alphaData)
end
if ~isempty(get(Img, 'XData')) && ~isempty(get(Img, 'YData'))
    set(IconAxes, ...
        'XLim', get(Img, 'XData')+[-0.5, 0.5], ...
        'YLim', get(Img, 'YData')+[-0.5, 0.5] ...
        );
end

set(IconAxes, ...
    'Visible', 'off', ...
    'YDir', 'reverse' ...
    );

movegui(d, 'center');
d.Visible = true;
end

%--------------------------------------------------------------------------
function buttonFitToText(thisBtn, minXpos)
% Check for text clipping and adjust position if needed

overhang = thisBtn.Extent(3) - thisBtn.Position(3);
if overhang >= 0
    % Split difference between x-position and width
    buffer = 5;
    thisBtn.Position(1) = max(thisBtn.Position(1)-floor(overhang/2), minXpos);
    overhang = thisBtn.Extent(3) - thisBtn.Position(3);
    thisBtn.Position(3) = thisBtn.Position(3) + ceil(overhang) + buffer;
end
end

%--------------------------------------------------------------------------
function openOptimizeTaskCallback(dlgHndl)
% Launch template Live Script

% Close dialog
delete(dlgHndl);

scgeatool;
end
