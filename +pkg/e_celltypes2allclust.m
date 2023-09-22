function [sce] = e_celltypes2allclust(sce, speciestag, mergec)

if nargin < 3, mergec = true; end
if nargin < 2, speciestag = 'mouse'; end

[c, cL] = grp2idx(sce.c_cluster_id);
organtag = "all";
databasetag = "panglaodb";
cLdisp = cL;
for i = 1:max(c)
    fprintf('Processing cluster....%d of %d\n', i, max(c));
    ptsSelected = c == i;
    [Tct] = pkg.local_celltypebrushed(sce.X, sce.g, ...
        sce.s, ptsSelected, ...
        speciestag, organtag, databasetag);
    if isempty(Tct)
        ctxt = 'Unknown';
    else
        ctxt = Tct.C1_Cell_Type{1};
    end
    ctxtdisp = strrep(ctxt, '_', '\_');
    ctxtdisp = sprintf('%s_{%d}', ctxtdisp, i);
    cLdisp{i} = ctxtdisp;

    ctxt = sprintf('%s_{%d}', ctxt, i);
    cL{i} = ctxt;
end
sce.c_cell_type_tx = string(cL(c));

if mergec
    [sce] = pkg.i_mergeSubCellNames(sce);
end
end