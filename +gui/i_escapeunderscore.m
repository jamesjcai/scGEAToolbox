function txt = i_escapeunderscore(txt)
%I_ESCAPEUNDERSCORE Escape underscores for TeX-interpreted graphics labels.
%   txt = I_ESCAPEUNDERSCORE(txt) replaces "_" with "\_" so titles, axis
%   tick labels, and legends render underscores literally instead of as TeX
%   subscripts.
%
%   A "_" that opens a complete subscript group - "_{...}" with a matching
%   closing brace - is left alone, because the toolbox writes cell subtypes
%   that way on purpose: "T cells_{Effector memory}" and "T cells_{1}" are
%   meant to render with the qualifier as a subscript. Every other
%   underscore, including one opening an unclosed "_{", is escaped.
%
%   A space is inserted before a kept subscript group, so the qualifier
%   does not crowd the word it hangs off: "T cells_{Effector memory}"
%   renders as "T cells" followed by a gap and then the subscript.
%
%   Unlike a bare STRREP call, this accepts any common label container:
%   char vectors, string arrays, cell arrays of char vectors, and cell
%   arrays of string scalars (which STRREP rejects with the error
%   "Cell elements must be character vectors"). Cell inputs return a cell
%   array of char vectors; all other inputs return a string array.

% A balanced, brace-delimited subscript group gets a space in front of it
% (unless one is already there), then every "_" outside such a group is
% escaped. Spacing first keeps the group recognisable to the second step.
% REGEXPREP will not substitute on a zero-width match, so the character
% before the group is captured and put back rather than looked behind.
gap = '(\S)(_\{[^{}]*\})';
plainunderscore = '_(?!\{[^{}]*\})';

wascell = iscell(txt);
txt = regexprep(string(txt), gap, '$1 $2');
txt = regexprep(txt, plainunderscore, '\\_');
if wascell
    txt = cellstr(txt);
end
end
