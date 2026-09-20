function s = i_plural(n, singular, plural)
%I_PLURAL A count and its noun, agreeing in number.
%
%   s = pkg.i_plural(n, singular)
%   s = pkg.i_plural(n, singular, plural)
%
%   pkg.i_plural(3, 'cell')                 "3 cells"
%   pkg.i_plural(1, 'cell')                 "1 cell"
%   pkg.i_plural(0, 'cell')                 "0 cells"
%   pkg.i_plural(2, 'analysis', 'analyses') "2 analyses"
%
% Returns a string. PLURAL defaults to SINGULAR with an "s"; pass it for a
% noun that does not pluralize that way.
%
% For counts written into a dialog the user reads. "1 cells" in a list of
% proposed merges is the sort of thing that makes a tool look unfinished, and
% "cell(s)" everywhere reads like a form rather than a sentence - worth the
% one call when the count can genuinely be one.
%
% see also: gui.callback_CollapseCellSubtypes, gui.callback_SubtypeAnnotation

singular = string(singular);
if nargin < 3, plural = singular + "s"; end
plural = string(plural);

s = string(n) + " " + plural;
s(n == 1) = string(n(n == 1)) + " " + singular;
end
