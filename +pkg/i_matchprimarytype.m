function [primary, ismatched] = i_matchprimarytype(labels, primarytypes)
% I_MATCHPRIMARYTYPE - Map annotation labels onto primary cell types.
%
%   [primary, ismatched] = pkg.i_matchprimarytype(labels, primarytypes)
%
% Inputs:
%   labels       - observed cell type labels, e.g. sce.c_cell_type_tx
%   primarytypes - the types a marker table knows about, e.g. the CellType
%                  column of cellsubtypes.xlsx
%
% Outputs:
%   primary   - for each label, the primary type it belongs to, "" when none
%   ismatched - primary ~= ""
%
% Gross annotations rarely name a type the way a marker table does. The labels
% SCE.ASSIGNCELLTYPE writes carry a cluster index ("T cells_{3}"), a lineage
% prefix ("CD8+ T cells"), a state ("T memory cells"), or a synonym
% ("T lymphocytes"), and an exact string comparison against "T cells" misses
% every one of them. This maps all of those onto the primary type so that
% subtype annotation can be offered for them.
%
% Matching is deliberately conservative - a label reaches a primary type only
% through one of these, in order:
%
%   1. the same name, once case, punctuation and a trailing cluster index are
%      normalized away
%   2. a known synonym of the primary type ("Tregs", "CTL", "T lymphocyte")
%   3. the primary type appearing as whole words inside the label, so
%      "CD8+ T cells" matches "T cells" while "Mast cells" does not
%   4. the label opening with the primary type's first word and closing with
%      its last, which is the "T memory cells" shape
%
% Nothing matches on a shared head noun alone: "Natural killer cells" is not
% "T cells", however PanglaoDB happens to file its subtypes.
%
% see also: pkg.i_normalizetypename, sc_csubtypeanno,
%           gui.callback_SubtypeAnnotation

% string() before (:), not after: 'T cells' as a char row vector would
% otherwise be reshaped into seven one-character types.
labels = string(labels);
primarytypes = string(primarytypes);
primarytypes = primarytypes(:);
primarytypes = primarytypes(strlength(primarytypes) > 0);

primary = strings(size(labels));
ismatched = false(size(labels));
if isempty(labels) || isempty(primarytypes), return; end

% One key per distinct label: annotations repeat over thousands of cells and
% the work below is per distinct name, not per cell.
[ulabels, ~, back] = unique(labels(:));
ukeys = in_key(ulabels);
pkeys = in_key(primarytypes);

uprimary = strings(size(ulabels));
for k = 1:numel(ukeys)
    if strlength(ukeys(k)) == 0, continue; end
    hit = in_findprimary(ukeys(k), pkeys);
    if hit > 0
        uprimary(k) = primarytypes(hit);
    end
end

primary(:) = uprimary(back);
ismatched = strlength(primary) > 0;
end

function idx = in_findprimary(labelkey, pkeys)
% First primary type the label reaches, 0 for none. The rules are tried in
% order of how specific they are, so a label that names one type outright is
% never claimed by another type's looser rule.

idx = 0;

hit = find(pkeys == labelkey, 1);
if ~isempty(hit), idx = hit; return; end

hit = find(arrayfun(@(p) ismember(labelkey, in_synonyms(p)), pkeys), 1);
if ~isempty(hit), idx = hit; return; end

hit = find(arrayfun(@(p) in_containswords(labelkey, p), pkeys), 1);
if ~isempty(hit), idx = hit; return; end

hit = find(arrayfun(@(p) in_bracketsword(labelkey, p), pkeys), 1);
if ~isempty(hit), idx = hit; return; end
end

function key = in_key(s)
% Normalized comparison key: lowercase, punctuation to spaces, the cluster
% index SCE.ASSIGNCELLTYPE appends dropped, and the last word singularized so
% "T cells" and "T cell" are one key.

key = pkg.i_normalizetypename(s);
key = regexprep(key, '\s+\d+$', '');            % 'T cells_{3}' -> 't cells'
key = regexprep(key, '(\w)s$', '$1');           % 't cells' -> 't cell'
key = regexprep(key, '\<lymphocyte$', 'cell');  % 'B lymphocytes' -> 'b cell'
key = strtrim(key);
end

function tf = in_containswords(labelkey, pkey)
% Does the primary type appear in the label as whole words? The word boundary
% is what keeps 'mast cell' out of 't cell'.

% Built with string concatenation, not [ ]: bracketing a string between two
% char literals builds a 1x3 string array of patterns rather than one pattern,
% and regexp then matches against any of them.
pat = "(^|\s)" + regexptranslate('escape', pkey) + "($|\s)";
tf = ~isempty(regexp(labelkey, pat, 'once'));
end

function tf = in_bracketsword(labelkey, pkey)
% Does the label open with the primary type's first word and close with its
% last, with something in between? That is the 'T memory cells' shape, which
% no amount of substring matching finds.

tf = false;
ptok = split(pkey, ' ');
if numel(ptok) < 2, return; end

ltok = split(labelkey, ' ');
if numel(ltok) <= numel(ptok), return; end

tf = ltok(1) == ptok(1) && ltok(end) == ptok(end);
end

function syn = in_synonyms(pkey)
% Synonyms of a primary type, as normalized keys. Only names that mean the
% type itself belong here: a subtype name is listed when it is commonly used
% as the whole population's label in a gross annotation (Tregs, CTLs), which
% is what a subtype run is meant to resolve.

persistent tbl
if isempty(tbl)
    tbl = containers.Map('KeyType', 'char', 'ValueType', 'any');
    tbl('t cell') = ["t", "t lymphocyte", "tcell", "tconv", ...
        "conventional t cell", "alpha beta t cell", "ab t cell", ...
        "cd3 t cell", "cd4 t cell", "cd8 t cell", "cd4 t", "cd8 t", ...
        "cd4 positive t cell", "cd8 positive t cell", "cd4", "cd8", ...
        "treg", "t reg", "ctl", "cytotoxic lymphocyte", "th1", "th2", ...
        "th17", "tfh", "tph", "tem", "tcm", "temra", "trm", "mait", ...
        "mait cell", "nkt", "nkt cell", "gamma delta", "gd t cell", ...
        "tnk", "t nk cell"];
    tbl('b cell') = ["b", "b lymphocyte", "bcell", "breg"];
    tbl('nk cell') = ["nk", "natural killer", "natural killer cell", ...
        "nk lymphocyte", "cd56 nk cell", "innate lymphoid cell 1"];
    tbl('monocyte') = ["mono", "cd14 mono", "cd16 mono", "monocytic cell"];
    tbl('dendritic cell') = ["dc", "dendritic", "cdc", "cdc1", "cdc2", ...
        "pdc", "modc", "mregdc", "myeloid dendritic cell"];
    tbl('neutrophil') = ["pmn", "neutrophil granulocyte", ...
        "polymorphonuclear cell", "polymorphonuclear leukocyte"];
    tbl('fibroblast') = ["caf", "cancer associated fibroblast", ...
        "fibroblastic cell", "fibroblastic stromal cell"];
    tbl('endothelial cell') = ["ec", "endothelium", "endothelial", ...
        "vascular endothelium"];
    tbl('epithelial cell') = ["epithelium", "epithelial"];
    tbl('astrocyte') = ["astroglia", "astroglial cell"];
    tbl('microglia') = ["microglial cell", "microglial", "microglium"];
    tbl('oligodendrocyte') = ["oligo", "oligodendroglia", ...
        "oligodendroglial cell", "opc", "oligodendrocyte precursor cell", ...
        "oligodendrocyte progenitor cell"];
    % Microglia are their own primary type in the marker table, so they are
    % deliberately not listed here; Kupffer cells and alveolar macrophages are
    % not, and have nowhere else to go.
    tbl('macrophage') = ["mphi", "mf", "macrophage cell", "kupffer cell", ...
        "alveolar macrophage", "monocyte derived macrophage", "tam", ...
        "tumor associated macrophage", "histiocyte"];
    tbl('neuron') = ["neurone", "neuronal cell", "nerve cell", ...
        "neuronal population", "excitatory neuron", "inhibitory neuron"];
end

syn = strings(0, 1);
if isKey(tbl, char(pkey))
    syn = tbl(char(pkey));
end
end
