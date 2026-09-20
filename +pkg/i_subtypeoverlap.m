function [primary, subtype, ismatched] = i_subtypeoverlap(labels)
%I_SUBTYPEOVERLAP Read a label that already names a cell subtype.
%
%   [primary, subtype, ismatched] = pkg.i_subtypeoverlap(labels)
%
% Inputs:
%   labels  - observed cell type labels, e.g. sce.c_cell_type_tx
%
% Outputs:
%   primary   - for each label, the primary type of PKG.I_SUBTYPEOVERLAPTABLE
%               it belongs to, "" when the label is not a subtype-level name
%   subtype   - the subtype that label already names, "" likewise
%   ismatched - primary ~= ""
%
% celltypes.xlsx and cellsubtypes.xlsx overlap: 37 names in the first are
% subtypes of the 14 primaries the second covers. A primary annotation run can
% therefore hand back "T regulatory cells" or "Plasma cells", which are not
% unrelated types but T cells and B cells whose subtype is already decided.
% This says so, by exactly the same key PKG.I_MATCHPRIMARYTYPE compares on, so
% "T regulatory cells", "T-regulatory cell" and "T regulatory cells_{2}" all
% read the same.
%
% Two callers want opposite things from the answer and both are right:
% PKG.I_MATCHPRIMARYTYPE uses PRIMARY to recognize the cell, and
% SC_CSUBTYPEANNO uses ISMATCHED to leave it alone - re-clustering a cell that
% already carries a subtype would only overwrite a finer label with a coarser
% guess.
%
% see also: pkg.i_subtypeoverlaptable, pkg.i_matchprimarytype,
%           sc_csubtypeanno, pkg.i_subtypecandidates

labels = string(labels);
primary = strings(size(labels));
subtype = strings(size(labels));
ismatched = false(size(labels));
if isempty(labels), return; end

T = pkg.i_subtypeoverlaptable();
tkeys = pkg.i_typekey(T.CellType);

% One lookup per distinct label: annotations repeat over thousands of cells.
[ulabels, ~, back] = unique(labels(:));
[tf, loc] = ismember(pkg.i_typekey(ulabels), tkeys);

uprimary = strings(size(ulabels));
usubtype = strings(size(ulabels));
uprimary(tf) = T.Primary(loc(tf));
usubtype(tf) = T.SubType(loc(tf));

primary(:) = uprimary(back);
subtype(:) = usubtype(back);
ismatched = strlength(primary) > 0;
end
