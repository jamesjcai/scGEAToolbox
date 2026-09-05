function c = i_accentred(fig)
% I_ACCENTRED  A red for node labels that stays readable on either theme.
%
%   c = gui.i_accentred(fig) returns the RGB triplet to use for the labels
%   a network view picks out in red - transcription factors in
%   gui.i_singlegraph, seed genes in gui.i_multigraphs.
%
%   These labels mark the nodes the reader opened the figure to find, so
%   they are the worst ones to leave low-contrast. A single hardcoded red
%   cannot serve both grounds: measured against the axes background the
%   viewers actually paint (white, and 0.071 grey in dark mode), pure red
%   gives a contrast ratio of 4.00:1 on light - under the 4.5:1 WCAG AA
%   floor - while a dark red deep enough to fix that falls to 2.87:1 on
%   dark. The two returned here clear the floor on their own ground:
%   6.52:1 and 7.08:1.
%
%   See also gui.i_getthemebkgcolor, gui.i_singlegraph, gui.i_multigraphs.

if nargin < 1 || isempty(fig), fig = gcf; end

if all(gui.i_getthemebkgcolor(fig) > 0.5)
    c = [0.75, 0, 0];
else
    c = [1, 0.45, 0.45];
end
end
