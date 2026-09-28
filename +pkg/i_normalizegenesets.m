function [setmatrx, setnames, setgenes] = i_normalizegenesets(setmatrx, setnames, setgenes)
%I_NORMALIZEGENESETS  Bring a gene set collection to the membership-matrix form.
%
%   [setmatrx, setnames, setgenes] = PKG.I_NORMALIZEGENESETS(setmatrx, ...
%   setnames, setgenes) accepts either of the two ways a collection is
%   passed around this toolbox and returns the first of them:
%
%     - the PKG.E_GETGENESETS triple: an nSets-by-nSetGenes membership
%       matrix, the set names, and the gene symbols labelling its columns;
%     - an nSets-by-1 cell array of gene-name lists, with SETGENES ignored
%       and built here from the union of those lists.
%
%   Missing SETNAMES are filled in as "Set1", "Set2", ... so that a caller
%   trying one ad-hoc set does not have to invent a name for it.
%
%   Shared by SC_GSETTEST and RUN.ML_GSEA so that the two agree on what they
%   accept. They test the same collections under different nulls, and a
%   collection one of them takes and the other refuses would be a trap.
%
% See also PKG.E_GETGENESETS, SC_GSETTEST, RUN.ML_GSEA.

if iscell(setmatrx) && ~isempty(setmatrx) && ~isnumeric(setmatrx{1})
    lists = cellfun(@(v) string(v(:)), setmatrx(:), UniformOutput=false);
    setgenes = unique(vertcat(lists{:}));
    setgenes = setgenes(strlength(setgenes) > 0);
    mat = false(numel(lists), numel(setgenes));
    for k = 1:numel(lists)
        mat(k, :) = ismember(setgenes, lists{k});
    end
    setmatrx = mat;
end
if isempty(setgenes)
    error("pkg:i_normalizegenesets:NoSetGenes", ...
        "SETGENES is empty. Pass the gene symbols for the columns of " + ...
        "SETMATRX, or pass SETMATRX as a cell array of gene-name lists.");
end
setgenes = string(setgenes(:));
if isempty(setnames)
    setnames = "Set" + string((1:size(setmatrx, 1))');
end
setnames = string(setnames(:));
if numel(setnames) ~= size(setmatrx, 1)
    error("pkg:i_normalizegenesets:SizeMismatch", ...
        "SETNAMES (%d) must have one entry per row of SETMATRX (%d).", ...
        numel(setnames), size(setmatrx, 1));
end
if numel(setgenes) ~= size(setmatrx, 2)
    error("pkg:i_normalizegenesets:SizeMismatch", ...
        "SETGENES (%d) must have one entry per column of SETMATRX (%d).", ...
        numel(setgenes), size(setmatrx, 2));
end

end
