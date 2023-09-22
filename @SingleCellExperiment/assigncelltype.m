function obj = assigncelltype(obj, speciesid)
% Assigning cell type identity to clusters.
if nargin < 2, speciesid = 'mouse'; end
[c, cL] = grp2idx(obj.c_cluster_id);
organtag = "all";
databasetag = "panglaodb";
for i = 1:max(c)
    ptsSelected = c == i;
    [Tct] = pkg.local_celltypebrushed(obj.X, obj.g, ...
        obj.s, ptsSelected, ...
        speciesid, organtag, databasetag);
    ctxt = Tct.C1_Cell_Type{1};
    ctxt = sprintf('%s_{%d}', ctxt, i);
    cL{i} = ctxt;
end
obj.c_cell_type_tx = string(cL(c));
end
