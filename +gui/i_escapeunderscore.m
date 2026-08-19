function txt = i_escapeunderscore(txt)
%I_ESCAPEUNDERSCORE Escape underscores for TeX-interpreted graphics labels.
%   txt = I_ESCAPEUNDERSCORE(txt) replaces every "_" with "\_" so titles,
%   axis tick labels, and legends render underscores literally instead of
%   as TeX subscripts.
%
%   Unlike a bare STRREP call, this accepts any common label container:
%   char vectors, string arrays, cell arrays of char vectors, and cell
%   arrays of string scalars (which STRREP rejects with the error
%   "Cell elements must be character vectors"). Cell inputs return a cell
%   array of char vectors; all other inputs return a string array.

if iscell(txt)
    txt = cellstr(strrep(string(txt), '_', '\_'));
else
    txt = strrep(string(txt), '_', '\_');
end
end
