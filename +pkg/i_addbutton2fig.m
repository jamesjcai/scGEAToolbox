function [pt] = i_addbutton2fig(toolbarHdl, sepTag, callbackFnc, imgFil, tooltipTxt)

mfolder = fileparts(mfilename('fullpath'));

if nargin < 2 || isempty(sepTag), sepTag = 'off'; end
if nargin < 3 || isempty(callbackFnc), callbackFnc = "version"; end
if nargin < 4 || isempty(imgFil), imgFil = 'list.gif'; end
if nargin < 5 || isempty(tooltipTxt), tooltipTxt = "Test"; end

if ischar(callbackFnc) || isstring(callbackFnc)
    callbackFnc = str2func(callbackFnc);
end

pt = uipushtool(toolbarHdl, 'Separator', sepTag);
try
    [img, map] = imread(fullfile(mfolder, '..', 'resources', 'Images', imgFil));
    ptImage = ind2rgb(img, map);
catch
    ptImage = rand(16, 16, 3);
end

pt.CData = ptImage;
pt.Tooltip = tooltipTxt;
pt.ClickedCallback = callbackFnc;
end