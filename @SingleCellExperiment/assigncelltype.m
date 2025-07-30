function obj = assigncelltype(obj, speciesid, keepclusterid)

% Assigning cell type identity to clusters.

if nargin < 3 || isempty(keepclusterid), keepclusterid = true; end
if nargin < 2 || isempty(speciesid), speciesid = 'human'; end

[c, cL] = findgroups(string(obj.c_cluster_id));
organtag = "all";
databasetag = "panglaodb";
for ik = 1:max(c)
    ptsSelected = c == ik;
    [Tct] = pkg.local_celltypebrushed(obj.X, obj.g, ...
        obj.s, ptsSelected, ...
        speciesid, organtag, databasetag);
    ctxt = Tct.C1_Cell_Type{1};
    if keepclusterid
        ctxt = sprintf('%s_{%d}', ctxt, ik);
    end
    cL(ik) = ctxt;
end
obj.c_cell_type_tx = cL(c);
end

