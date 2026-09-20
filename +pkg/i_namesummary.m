function [s] = i_namesummary(names, maxshown)
%I_NAMESUMMARY A list of names short enough to put in a dialog.
%
%   s = pkg.i_namesummary(names)
%   s = pkg.i_namesummary(names, maxshown)
%
%   "T cells", "T cells and B cells", "T cells, B cells and 2 others" -
%   enough to recognise what is being talked about without a dialog that
%   scrolls. MAXSHOWN is how many are named before the count takes over
%   (default 2).
%
%   See also gui.i_getcustommarkers, gui.callback_ExpandBrushedCells.

if nargin < 2 || isempty(maxshown), maxshown = 2; end

names = strtrim(string(names));
names = names(strlength(names) > 0);
names = names(:).';

if isempty(names)
    s = "";
elseif isscalar(names)
    s = names;
elseif numel(names) <= maxshown + 1
    % One more than MAXSHOWN still reads better named than counted:
    % "A, B and C" beats "A, B and 1 other".
    s = strjoin(names(1:end-1), ", ") + " and " + names(end);
else
    s = strjoin(names(1:maxshown), ", ") + ...
        sprintf(" and %d others", numel(names)-maxshown);
end
end
