function thisc = i_crossgroups(parts, sep, missingtag)
%I_CROSSGROUPS  Cross several per-cell grouping vectors into one label.
%
%   thisc = gui.i_crossgroups(parts)
%   thisc = gui.i_crossgroups(parts, sep, missingtag)
%
% PARTS is either an n-by-k array (any type STRING accepts) or a 1-by-k cell
% array of n-element vectors - one grouping variable per column/element.
% Returns an n-by-1 string in which each cell carries every variable's
% level joined by SEP, e.g. "Macrophages | IL".
%
% Only combinations that actually occur end up in the result: crossing is
% per cell, not a product of levels, so a combination nothing falls into
% never becomes an empty group downstream.
%
% SEP defaults to " | ". It avoids the TeX-special characters (_ ^ { } \)
% that the plotting layer escapes, so composite labels survive being used
% as axis ticks.
%
% MISSINGTAG (default "(n/a)") replaces missing and empty levels BEFORE the
% join. This is not cosmetic. A <missing> element propagates through JOIN
% and turns the entire composite for that cell into <missing>, which
% FINDGROUPS then assigns to a NaN group that GRPSTATS and CATEGORICAL
% quietly drop - so cells would vanish from a plot rather than show up as
% an identifiable group. An empty string does not propagate but renders as
% an invisible gap, which is worse for being unnoticeable.
%
% See also GUI.I_SELECTNCLASS, GUI.I_SELECT1CLASS.

if nargin < 3 || isempty(missingtag), missingtag = "(n/a)"; end
if nargin < 2 || isempty(sep), sep = " | "; end

if iscell(parts)
    k = numel(parts);
    assert(k > 0, 'gui:i_crossgroups:empty', 'No grouping variable given.');
    n = numel(parts{1});
    S = strings(n, k);
    for j = 1:k
        v = parts{j};
        assert(numel(v) == n, 'gui:i_crossgroups:length', ...
            'Grouping variable %d has %d elements, expected %d.', ...
            j, numel(v), n);
        S(:, j) = string(v(:));
    end
else
    S = string(parts);
    assert(~isempty(S), 'gui:i_crossgroups:empty', 'No grouping variable given.');
end

S(ismissing(S) | strlength(S) == 0) = missingtag;
thisc = join(S, sep, 2);
end
