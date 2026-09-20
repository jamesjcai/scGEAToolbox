function T = i_subtypeoverlaptable()
%I_SUBTYPEOVERLAPTABLE Names that assets/PanglaoDB's two tables both cover.
%
%   T = pkg.i_subtypeoverlaptable()
%
% Returns a table with one row per cell type name that appears in
% assets/PanglaoDB/celltypes.xlsx and names, in different words, a row of
% assets/PanglaoDB/cellsubtypes.xlsx:
%
%   CellType  - the name as celltypes.xlsx spells it, e.g. "T regulatory cells"
%   Primary   - the CellType column of cellsubtypes.xlsx, e.g. "T cells"
%   SubType   - the SubType column of cellsubtypes.xlsx, e.g. "Regulatory"
%
% The two workbooks were curated apart and overlap: celltypes.xlsx is meant to
% hold primary types and carries 37 names that are really subtypes of the 14
% primaries cellsubtypes.xlsx covers. PKG.E_DETERMINECELLTYPE scores a cluster
% against all of celltypes.xlsx at once, so one primary annotation run can
% return "T cells" for one cluster and "T regulatory cells" for the next. This
% table is what lets the rest of the toolbox read that second label for what it
% is - a T cell that already has its subtype - instead of an unrelated type.
%
% It is deliberately a code-side lookup: neither workbook is modified, so the
% markers behind both spellings stay available and a primary annotation run
% still returns exactly what it returned before.
%
% Only names whose subtype the bundled table actually lists are here. Cell
% types that are biologically a subdivision of a covered primary but have no
% row to map onto - "Thymocytes", "Pyramidal cells",
% "Myeloid-derived suppressor cells" - are left out, because there is nothing
% to say about them beyond which primary they belong to, and
% PKG.I_MATCHPRIMARYTYPE already answers that.
%
% A spelling that celltypes.xlsx has since corrected stays here rather than
% being dropped with it. Labels outlive the workbook: sce.c_cell_type_tx is
% saved into .mat files, so a dataset annotated by an earlier release still
% carries the old name and still has to be readable. Such a row is marked
% below, and TESTS/SUBTYPEOVERLAPTEST lists the same names so that the check
% tying this table to the workbook can tell a retired spelling from a typo.
%
% see also: pkg.i_subtypeoverlap, pkg.i_matchprimarytype, sc_csubtypeanno

persistent cached
if ~isempty(cached)
    T = cached;
    return;
end

% CellType (as celltypes.xlsx spells it), Primary, SubType.
rows = [ ...
    "Adrenergic neurons",               "Neurons",           "Adrenergic"
    "Cholinergic neurons",              "Neurons",           "Cholinergic"
    "Dopaminergic neurons",             "Neurons",           "Dopaminergic"
    "Enteric neurons",                  "Neurons",           "Enteric"
    "GABAergic neurons",                "Neurons",           "GABAergic"
    "Glutamatergic neurons",            "Neurons",           "Glutamatergic"
    "Glutaminergic neurons",            "Neurons",           "Glutamatergic"  % retired spelling
    "Glycinergic neurons",              "Neurons",           "Glycinergic"
    "Immature neurons",                 "Neurons",           "Immature"
    "Interneurons",                     "Neurons",           "Interneurons"
    "Motor neurons",                    "Neurons",           "Motor"
    "Noradrenergic neurons",            "Neurons",           "Noradrenergic"
    "Purkinje neurons",                 "Neurons",           "Purkinje"
    "Serotonergic neurons",             "Neurons",           "Serotonergic"
    "Trigeminal neurons",               "Neurons",           "Trigeminal"
    "T cells naive",                    "T cells",           "Naive"
    "T cytotoxic cells",                "T cells",           "Cytotoxic"
    "T follicular helper cells",        "T cells",           "Follicular helper"
    "T helper cells",                   "T cells",           "Helper"
    "T memory cells",                   "T cells",           "Memory"
    "T regulatory cells",               "T cells",           "Regulatory"
    "Natural killer T cells",           "T cells",           "Natural killer"
    "B cells naive",                    "B cells",           "Naive"
    "B cells memory",                   "B cells",           "Memory"
    "Plasma cells",                     "B cells",           "Plasma cell"
    "Plasmacytoid dendritic cells",     "Dendritic cells",   "Plasmacytoid"
    "Langerhans cells",                 "Dendritic cells",   "Langerhans"
    "Alveolar macrophages",             "Macrophages",       "Tissue-resident"
    "Kupffer cells",                    "Macrophages",       "Tissue-resident"
    "Red pulp macrophages",             "Macrophages",       "Tissue-resident"
    "Myofibroblasts",                   "Fibroblasts",       "Myofibroblastic"
    "Oligodendrocyte progenitor cells", "Oligodendrocytes",  "Precursor"
    "Basal cells",                      "Epithelial cells",  "Basal"
    "Ciliated cells",                   "Epithelial cells",  "Ciliated"
    "Goblet cells",                     "Epithelial cells",  "Goblet"
    "Airway goblet cells",              "Epithelial cells",  "Goblet"
    "Luminal epithelial cells",         "Epithelial cells",  "Luminal secretory"];

T = table(rows(:, 1), rows(:, 2), rows(:, 3), ...
    'VariableNames', {'CellType', 'Primary', 'SubType'});
cached = T;
end
