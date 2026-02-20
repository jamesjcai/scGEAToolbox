function [a] = i_setautumncolor(c, cmapname, rev, grayz, ax, parentfig)
% i_setautumncolor - Apply autumn colormap with auto dark-mode adjustments
%
% Inputs:
%   c        - data values (used to detect constant values)
%   cmapname - base colormap (default 'autumn')
%   rev      - reverse colormap (default true)
%   grayz    - set lowest color to gray (default true)
%   ax       - target axis (default gca)
%
% Output:
%   a        - final colormap matrix

if nargin < 6, parentfig = []; end
if nargin < 5, ax = []; end
if nargin < 4, grayz = true; end
if nargin < 3, rev = true; end
if nargin < 2, cmapname = 'autumn'; end
if isempty(ax), ax = gca; end

% isDark = false;
if isempty(parentfig)
    % Detect dark mode from figure background
    figBG = get(groot, 'DefaultFigureColor');
    isDark = mean(figBG) < 0.5;   % heuristic: dark mode if average brightness < 0.5
else
    try
        if strcmp('light', parentfig.Theme.BaseColorStyle)
            isDark = false;
        else
            isDark = true;  % dark mode detected
        end
    catch
        isDark = false;  % default to light mode if error occurs
    end    
end

% start from base colormap
a = feval(cmapname, 256);

% reverse if requested
if rev
    a = flipud(a);
end

% adjust autumn if in dark mode
if strcmpi(cmapname, 'autumn') && isDark
    % darken the bright yellow end toward gold
    a(:,2) = a(:,2) * 0.85;    % reduce green channel
    a(:,3) = a(:,3) * 0.6;     % reduce blue channel
end

% set gray for background/constant values
if grayz
    if isDark
        gcol = [0.5 0.5 0.5];  % medium gray for dark mode
    else
        gcol = [0.8 0.8 0.8];  % lighter gray for light mode
    end
    a(1,:) = gcol;
    if isscalar(unique(c))     % constant data → all gray
        a(:) = gcol(1);
    end
end

colormap(ax, a);
end

%{
function [a] = i_setautumncolor(c, cmapname, rev, grayz, ax)
if nargin < 5, ax = []; end
if nargin < 4, grayz = true; end
if nargin < 3, rev = true; end
if nargin < 2, cmapname = 'autumn'; end

if isempty(ax), ax = gca; end

colormap(ax, "default");
a = colormap(ax, cmapname);
if strcmpi(cmapname, 'autumn')
    if rev
        a = flipud(a);
    end
end
if grayz
    a(1, :) = [.8, .8, .8];
    if isscalar(unique(c)), a(:) = 0.8; end
end
colormap(ax, a);
end
%}

