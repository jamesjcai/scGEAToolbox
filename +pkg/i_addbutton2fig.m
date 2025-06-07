function [pt, ptImage] = i_addbutton2fig(toolbarHdl, sepTag, callbackFnc, imgFil, tooltipTxt)

% i_addbutton2fig - Add a button to a figure toolbar
%
% Syntax:
%   pt = i_addbutton2fig(toolbarHdl, sepTag, callbackFnc, imgFil, tooltipTxt)
%
% Inputs:
%   toolbarHdl - Handle to the target toolbar
%   sepTag - (Optional) 'on' or 'off' for visual separation (default: 'off')
%   callbackFnc - (Optional) Callback function or name (default: "version")
%   imgFil - (Optional) Image filename for button icon (default: empty)
%   tooltipTxt - (Optional) Tooltip text for the button (default: "Test")
%
% Outputs:
%   pt - Handle to the created push tool

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

        bgColor = get(toolbarHdl.Parent, 'Color');

        % bgColor = [0.9569    0.9569    0.9569];
        % bgColor = [0.9400    0.9400    0.9400];
        % transparentIcon = zeros([iconSize, 3]); % RGB values all set to zero (black)
        ptImage = repmat(reshape(bgColor, [1, 1, 3]), iconSize); % Match background
        pt = uipushtool(toolbarHdl, ...
            'CData', ptImage, ...
            'TooltipString', '');
    else
%        pt = uipushtool(toolbarHdl, 'Separator', sepTag);

        if ischar(imgFil) || isstring(imgFil)
            imgPath = fullfile(mfolder, '..', 'assets', 'Images', imgFil);
            if ~isfile(imgPath)
                % warning('Image file "%s" not found. Using a random icon.', imgPath);
            end
            try
                [ptImage, map] = imread(imgPath);
                if ~isempty(map), ptImage = ind2rgb(ptImage, map); end
            catch ME
                warning('Failed to load image "%s". Using random icon instead.\nError: %s', imgPath, ME.message);
                ptImage = rand(16, 16, 3);
            end
            if size(ptImage, 3) == 1, ptImage = cat(3, ptImage, ptImage, ptImage); end
        else
            ptImage = imgFil;
        end

        pt = uipushtool(toolbarHdl, ...
            'CData', ptImage, ...
            'TooltipString', tooltipTxt, ...
            'Separator', sepTag, 'ClickedCallback', callbackFnc);

        %        pt.CData = ptImage;
        %        pt.Tooltip = tooltipTxt;
        %        pt.ClickedCallback = callbackFnc;
    end
end


function resizedImage = resizeTo16x16(ptImage)
    % Resize an input image to 16x16 using interpolation
    % Input: 
    %   ptImage - RGB image (MxNx3) or grayscale image (MxN)
    % Output: 
    %   resizedImage - Resized image (16x16x3 or 16x16)
    
    % Get the size of the input image
    [rows, cols, channels] = size(ptImage);
    
    % Create the original grid
    [x, y] = meshgrid(1:cols, 1:rows);
    
    % Create the target grid for 16x16 resizing
    [xq, yq] = meshgrid(linspace(1, cols, 16), linspace(1, rows, 16));
    
    % Initialize the resized image
    if channels == 1
        % Grayscale image
        resizedImage = interp2(x, y, ptImage, xq, yq, 'linear');
    else
        % RGB image
        resizedImage = zeros(16, 16, channels);
        for channel = 1:channels
            resizedImage(:, :, channel) = interp2(x, y, ptImage(:, :, channel), xq, yq, 'linear');
        end
    end
end
