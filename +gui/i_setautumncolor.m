function [a] = i_setautumncolor(c, cmapname, rev, grayz)
if nargin < 4, grayz = true; end
if nargin < 3, rev = true; end
if nargin < 2, cmapname = 'autumn'; end
colormap default
a = colormap(cmapname);
if strcmpi(cmapname, 'autumn')
    if rev
        a = flipud(a);
    end
end
if grayz
    a(1, :) = [.8, .8, .8];
    if isscalar(unique(c)), a(:) = 0.8; end
end
colormap(a);
end