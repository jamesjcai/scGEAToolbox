function [parent, form] = i_labelparents(labels)
%I_LABELPARENTS The cell type each subtype label collapses back to.
%
%   [parent, form] = gui.i_labelparents(labels)
%
% Outputs, one element per label:
%   parent - the cell type it collapses to, or the label itself when it does
%            not collapse to anything
%   form   - "brace", "paren", "name", or "" when PARENT is the label
%
% Three ways a label can name a subtype, tried in that order:
%
%   brace  'T cells_{Regulatory}', 'T cells_{3}'. What this toolbox writes,
%          and nothing else does.
%   paren  'T cells (Regulatory)'. Also this toolbox - but celltypes.xlsx has
%          names of that shape too, 'Endothelial cells (aorta)', so a caller
%          has to confirm rather than apply it.
%   name   'Plasma cells', 'Interneurons', 'CD8+ T cells'. No suffix at all.
%
% The third is why this exists. Collapsing used to be suffix stripping and
% nothing else, so it could only undo labels this toolbox had written. A
% dataset annotated straight out of celltypes.xlsx carries no suffix
% anywhere, and celltypes.xlsx names 37 subtypes of the primaries
% cellsubtypes.xlsx covers - 'Plasma cells' is a B cell and 'Interneurons' is
% a neuron, with nothing in either string to say so. Those labels had no way
% back to their major type at all, which also left SC_CSUBTYPEANNO telling
% people to collapse a label that could not be collapsed.
%
% PKG.I_MATCHPRIMARYTYPE answers that, and answers more besides: it reaches a
% primary from a lineage prefix ('CD8+ T cells') or a synonym ('Tregs') as
% well as from the curated equivalences in PKG.I_SUBTYPEOVERLAP. Those are
% guesses from the shape of a name where the suffix forms are records of what
% was done, so a caller should offer them and let the user decide - which is
% what GUI.CALLBACK_COLLAPSECELLSUBTYPES does by leaving them unticked.
%
% see also: pkg.i_splitsubtypelabel, pkg.i_matchprimarytype,
%           pkg.i_subtypeoverlap, gui.callback_CollapseCellSubtypes

labels = string(labels);
[parent, ~, form] = pkg.i_splitsubtypelabel(labels);
if isempty(labels), return; end

% Only labels that gave up no suffix. One that did has already said what its
% parent is, and said it as a record rather than as a guess.
plain = parent == labels;
if ~any(plain), return; end

[~, ~, primarytypes] = pkg.i_subtypecandidates();
if isempty(primarytypes), return; end

primary = pkg.i_matchprimarytype(labels(plain), primarytypes);

% A label that IS a primary type reaches itself, which is not a collapse.
hit = strlength(primary) > 0 & primary ~= labels(plain);
if ~any(hit), return; end

at = find(plain);
parent(at(hit)) = primary(hit);
form(at(hit)) = "name";
end
