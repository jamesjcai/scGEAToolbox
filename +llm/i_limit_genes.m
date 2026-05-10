function sce = i_limit_genes(sce, n, mustkeep)
% I_LIMIT_GENES  Remove noise genes then keep top n HVGs.
%
%   sce = llm.i_limit_genes(sce, n)
%   sce = llm.i_limit_genes(sce, n, mustkeep)
%
%   Applies the same gene-reduction pipeline as gui.callback_Select5000Genes:
%     1. Remove mitochondrial genes     (sce.rmmtgenes)
%     2. Remove hemoglobin genes        (sce.rmhemoglobingenes)
%     3. Remove ribosomal genes         (sce.rmribosomalgenes)
%     4. Remove genes containing 'orf' or '-AS' (pseudogenes, antisense)
%     5. Remove genes starting with 'LINC' (lncRNAs)
%     6. Remove mouse predicted genes   (Gm[0-9]{3,})
%     7. Remove mouse RIKEN clones      (ends with 'Rik')
%     8. Keep top n HVGs via sc_splinefit
%
%   mustkeep - gene name that must survive all filters (e.g. KO gene).

if nargin < 3, mustkeep = ""; end

sce = sce.rmmtgenes;
sce = sce.rmhemoglobingenes;
sce = sce.rmribosomalgenes;

g = sce.g;
rm = contains(g, 'orf') | contains(g, '-AS') | contains(g, '-as') | ...
     startsWith(g, 'LINC') | ...
     ~cellfun(@isempty, regexp(g, 'Gm[0-9][0-9][0-9]')) | ...
     endsWith(g, 'Rik');
if strlength(mustkeep) > 0
    rm(g == mustkeep) = false;
end
if any(rm)
    sce.X(rm, :) = [];
    sce.g(rm) = [];
end

T_hvg = sc_splinefit(sce.X, sce.g);
glist = T_hvg.genes(1:min(n, sce.NumGenes));
if strlength(mustkeep) > 0 && ~ismember(mustkeep, glist)
    glist(end) = mustkeep;
end
[~, hidx] = ismember(glist, sce.g);
hidx = hidx(hidx > 0);
sce.X = sce.X(hidx, :);
sce.g = sce.g(hidx);
end
