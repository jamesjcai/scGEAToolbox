function obj = assigncelltype(obj, speciesid, keepclusterid, keepold)
%ASSIGNCELLTYPE Assign cell type identity to clusters from marker genes.
%
%   obj = obj.assigncelltype()
%   obj = obj.assigncelltype(speciesid, keepclusterid, keepold)
%
%   Matches each cluster in c_cluster_id against PanglaoDB and writes the
%   result to c_cell_type_tx.
%
%   speciesid      'human' (default) or 'mouse'.
%   keepclusterid  append '_{k}' to each label. Default true.
%   keepold        stash any existing cell type labels in a new numbered
%                  'old_cell_type_N' cell attribute before overwriting, so
%                  an earlier annotation is not silently lost. Default
%                  true. There is no dialog to ask on this path, so the
%                  flag is how a caller that means to discard them says so.
%
%   See also SC_ANNOTATECELLS, PKG.I_STASHCELLTYPEHISTORY.

if nargin < 4 || isempty(keepold), keepold = true; end
if nargin < 3 || isempty(keepclusterid), keepclusterid = true; end
if nargin < 2 || isempty(speciesid), speciesid = 'human'; end

[c, cL] = findgroups(string(obj.c_cluster_id));
organtag = "all";
databasetag = "panglaodb";
for ik = 1:max(c)
    ptsSelected = c == ik;
    [Tct] = pkg.i_celltypebrushed(obj.X, obj.g, ...
        obj.s, ptsSelected, ...
        speciesid, organtag, databasetag);
    ctxt = Tct.C1_Cell_Type{1};
    if keepclusterid
        ctxt = sprintf('%s_{%d}', ctxt, ik);
    end
    cL(ik) = ctxt;
end
if keepold
    pkg.i_stashcelltypehistory(obj);
end
obj.c_cell_type_tx = cL(c);
end
