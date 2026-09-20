function [missing, total] = i_missingmarkers(sce, Tm)
%I_MISSINGMARKERS Marker genes of a list that are not in the data.
%
%   [missing, total] = pkg.i_missingmarkers(sce, Tm)
%
%   Tm is a two-column marker table - the first column the type name, the
%   second its markers as one comma separated string - as
%   PKG.I_PARSEMARKERLIST returns. MISSING is the symbols that do not appear in
%   SCE.G, and TOTAL how many distinct symbols the list names.
%
%   Markers that are not in this dataset at all score nothing, so a type whose
%   whole list is missing can never be assigned. A typo or the wrong species is
%   the usual reason, and neither shows up anywhere else.
%
%   See also gui.i_warnmissingmarkers, pkg.i_parsemarkerlist.

missing = strings(0, 1);
total = 0;
if nargin < 2 || isempty(Tm) || isempty(sce), return; end

markers = strtrim(split(strjoin(upper(string(Tm{:, 2})), ","), ","));
markers = unique(markers(strlength(markers) > 0));
total = numel(markers);

missing = markers(~ismember(markers, upper(string(sce.g))));
end
