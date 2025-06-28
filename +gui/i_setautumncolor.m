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