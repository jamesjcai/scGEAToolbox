function label = i_embeddinglabel(sce)
%I_EMBEDDINGLABEL Name of the method behind the embedding SCE.S, if known.
%
%   label = pkg.i_embeddinglabel(sce) compares SCE.S with the embeddings
%   stored in SCE.STRUCT_CELL_EMBEDDINGS and returns the method of the one
%   it matches, as an axis label: "tSNE", "UMAP", "PHATE", "MetaViz" or
%   "Monocle". It returns "" when SCE.S matches none of them, or matches
%   one whose method this function does not know - for example after
%   Harmony has replaced SCE.S - so a caller can fall back to asking.
%
%   EMBEDCELLS stores SINGLE(SCE.S) while SCE.S stays double, so the
%   comparison is made in single precision.
%
%   See also SINGLECELLEXPERIMENT/EMBEDCELLS, PKG.E_MAKEEMBEDSTRUCT.

label = "";
embeddings = sce.struct_cell_embeddings;
if ~isstruct(embeddings) || isempty(sce.s), return; end

prefixes = ["tsne", "umap", "phate", "metaviz", "monocle"];
labels = ["tSNE", "UMAP", "PHATE", "MetaViz", "Monocle"];

s = single(sce.s);
names = string(fieldnames(embeddings));
for k = 1:numel(names)
    e = embeddings.(names(k));
    if isempty(e) || ~isequal(size(e), size(s)) || ~isequal(single(e), s)
        continue;
    end
    % One test per prefix: STARTSWITH with a pattern array answers only
    % whether ANY pattern matches.
    hit = arrayfun(@(p) startsWith(lower(names(k)), p), prefixes);
    if any(hit)
        label = labels(find(hit, 1));
        return;
    end
end
end
