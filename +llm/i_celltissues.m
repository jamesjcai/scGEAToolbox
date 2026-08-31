function tissues = i_celltissues()
%I_CELLTISSUES Tissues covered by the cell-type reference dictionary.
%
%   tissues = LLM.I_CELLTISSUES()
%
%   The names are hard-coded so a chooser can be shown before anything is
%   downloaded; LLM.I_CELLCONTEXT fetches the markers themselves on first
%   use. Any other tissue name still works as free text - it just does not
%   come with a reference cell-type list.
%
%   The CellTypeAI README advertises 16 tissues, but its data file carries a
%   17th, Blood, with 30 cell types. It is listed here because it is real,
%   well populated, and the most likely choice for PBMC data.
%
%   See also LLM.I_CELLCONTEXT, LLM.E_CELLTYPEANNO.

tissues = ["Adrenal", "Blood", "Brain", "Eye", "Heart", "Immune system", ...
    "Intestine", "Kidney", "Liver", "Lung", "Muscle", "Pancreas", ...
    "Placenta", "Skin", "Spleen", "Stomach", "Thymus"];
end
