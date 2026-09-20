function [wvalu, wgene, celltypev, markergenev] = i_markerweights(Tm)
%I_MARKERWEIGHTS Weights and gene lists for scoring against a marker table.
%
%   [wvalu, wgene, celltypev, markergenev] = pkg.i_markerweights(Tm)
%
%   Tm is a two-column marker table - the first column the type name, the
%   second its markers as one comma separated string - as
%   PKG.I_PARSEMARKERLIST returns. The four outputs are the arguments
%   PKG.E_DETERMINECELLTYPE takes after the cell selection.
%
%   Everything is upper-cased on the way out. PKG.E_DETERMINECELLTYPE splits
%   MARKERGENEV and compares each symbol with == against WGENE and against
%   upper(sce.g), so a marker arriving in mixed case matches neither and is
%   silently dropped from its type's score.
%
%   See also pkg.e_markerweight, pkg.e_determinecelltype,
%   pkg.i_parsemarkerlist, gui.i_getcustommarkers.

celltypev = string(Tm{:, 1});
markergenev = upper(string(Tm{:, 2}));

Tw = pkg.e_markerweight(table(celltypev, markergenev, ...
    'VariableNames', {'CellType', 'PositiveMarkers'}));

wvalu = Tw.Var2;
wgene = string(Tw.Var1);
end
