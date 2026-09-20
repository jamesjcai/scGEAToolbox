function [parent, subtype, form] = i_splitsubtypelabel(labels)
%I_SPLITSUBTYPELABEL Split a cell type label into its parent type and subtype.
%
%   [parent, subtype, form] = pkg.i_splitsubtypelabel(labels)
%
%   Recognizes the two suffixes this toolbox writes onto a cell type label:
%
%     'T cells_{3}'           SCE.ASSIGNCELLTYPE's cluster index
%     'T cells_{Regulatory}'  SC_CSUBTYPEANNO with formatid 2
%     'T cells (Regulatory)'  SC_CSUBTYPEANNO with formatid 1
%
%   Outputs, one element per input label:
%     parent   the label with its last suffix removed, or the label itself
%              when it carries none
%     subtype  what was removed, "" when nothing was
%     form     "brace", "paren", or "" - which suffix was found
%
%   Only the LAST suffix is removed, so 'T cells_{3} (Regulatory)' gives a
%   parent of 'T cells_{3}'. Peeling the next level is another call; doing it
%   in one pass would merge two different groupings in a single step without
%   the caller ever seeing the intermediate.
%
%   THE PARENTHESIZED FORM IS AMBIGUOUS AND THIS FUNCTION DOES NOT RESOLVE IT.
%   'Endothelial cells (aorta)' is a cell type name in celltypes.xlsx, not a
%   subtype of 'Endothelial cells'. The collision is not hypothetical:
%   cellsubtypes.xlsx lists Arterial, Venous, Capillary, Lymphatic and Tip
%   cell for that same primary, so 'Endothelial cells (Arterial)' - which
%   SC_CSUBTYPEANNO really does write - has exactly the shape of the name
%   above, and nothing in either string says which of the two it is. Callers
%   must show the split and let the user confirm it rather than applying it
%   blind - see GUI.CALLBACK_COLLAPSECELLSUBTYPES, which preselects only the
%   brace form for exactly this reason.
%
%   See also sc_csubtypeanno, gui.callback_CollapseCellSubtypes,
%   SingleCellExperiment/assigncelltype.

labels = string(labels);
parent = labels;
subtype = strings(size(labels));
form = strings(size(labels));
if isempty(labels), return; end

% One match per distinct label: annotations repeat over thousands of cells.
[ulabels, ~, back] = unique(labels(:));
uparent = ulabels;
usubtype = strings(size(ulabels));
uform = strings(size(ulabels));

for k = 1:numel(ulabels)
    [uparent(k), usubtype(k), uform(k)] = in_split(ulabels(k));
end

parent(:) = uparent(back);
subtype(:) = usubtype(back);
form(:) = uform(back);
end

function [p, s, f] = in_split(label)
p = label;
s = "";
f = "";

t = strtrim(label);
if strlength(t) == 0, return; end

% Greedy on the parent, so the LAST group is the one taken. Anchored at the
% end, so a group in the middle of a name is left alone.
tok = regexp(t, '^(.*\S)_\{(.+)\}$', 'tokens', 'once');
if ~isempty(tok)
    f = "brace";
else
    tok = regexp(t, '^(.*\S)\s*\((.+)\)$', 'tokens', 'once');
    if ~isempty(tok)
        f = "paren";
    end
end
if isempty(tok), return; end

candidate = strtrim(tok{1});
if strlength(candidate) == 0
    % '(Regulatory)' on its own: there is no parent to fall back to, so this
    % is a cell type whose name happens to be parenthesized.
    f = "";
    return;
end

p = candidate;
s = strtrim(tok{2});
end
