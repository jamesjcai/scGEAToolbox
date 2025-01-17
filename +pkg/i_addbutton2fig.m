function [pt] = i_addbutton2fig(toolbarHdl, sepTag, callbackFnc, imgFil, tooltipTxt)

mfolder = fileparts(mfilename('fullpath'));

if nargin < 2 || isempty(sepTag), sepTag = 'off'; end
if nargin < 3 || isempty(callbackFnc), callbackFnc = "version"; end
if nargin < 4, imgFil = []; end
if nargin < 5 || isempty(tooltipTxt), tooltipTxt = "Test"; end

if ischar(callbackFnc) || isstring(callbackFnc)
    callbackFnc = str2func(callbackFnc);
end

    if isempty(imgFil)
        % % Add a "spacer" by using an invisible push tool
        iconSize = [16, 16]; 

        %{
        blankIcon = zeros([iconSize, 3]); % Black icon (invisible on dark themes)
        pt = uipushtool(toolbarHdl, ...
            'CData', blankIcon, ...
            'TooltipString', '');
        % pt = uipushtool(toolbarHdl, ...
        %     'Separator', 'off', ...  % Adds visual separation
        %     'Visible', 'on');      % Hides this tool to act as a spacer  
        %}

        
        % bgColor = get(hFig, 'Color'); % Figure background color matches toolbar

        bgColor = [0.9400    0.9400    0.9400];
        % transparentIcon = zeros([iconSize, 3]); % RGB values all set to zero (black)
        transparentIcon = repmat(reshape(bgColor, [1, 1, 3]), iconSize); % Match background
        pt = uipushtool(toolbarHdl, ...
            'CData', transparentIcon, ...
            'TooltipString', '');
    else
        pt = uipushtool(toolbarHdl, 'Separator', sepTag);
        try
            [ptImage, map] = imread(fullfile(mfolder, '..', 'resources', 'Images', imgFil));
            if ~isempty(map)
                ptImage = ind2rgb(ptImage, map);
            end
        catch
            ptImage = rand(16, 16, 3);
        end
        
        if size(ptImage, 3) == 1
        %     colormapName = gray(256);
        %     ptImage = ind2rgb(uint8(ptImage * 255), colormapName);
        ptImage = cat(3, ptImage, ptImage, ptImage);
        end
        
        pt.CData = ptImage;
        pt.Tooltip = tooltipTxt;
        pt.ClickedCallback = callbackFnc;
    end
end