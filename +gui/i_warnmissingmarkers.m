function [missing] = i_warnmissingmarkers(parentfig, sce, Tm, dlgtitle)
%I_WARNMISSINGMARKERS Warn about marker genes that are not in the data.
%
%   missing = gui.i_warnmissingmarkers(parentfig, sce, Tm)
%   missing = gui.i_warnmissingmarkers(parentfig, sce, Tm, dlgtitle)
%
%   Tm is a two-column marker table - Var1 the type name, Var2 its markers as
%   one comma separated string - as PKG.I_PARSEMARKERLIST returns. MISSING is
%   the marker symbols that do not appear in SCE.G, and a dialog names them.
%   PKG.I_MISSINGMARKERS does the finding; this adds the dialog.
%
%   Markers that are not in this dataset at all score nothing, so a type whose
%   whole list is missing can never be assigned. A typo or the wrong species is
%   the usual reason, and neither shows up anywhere else.
%
%   See also pkg.i_missingmarkers, gui.i_getcustommarkers,
%   gui.callback_SubtypeAnnotation.

if nargin < 4 || isempty(dlgtitle), dlgtitle = 'Customized Marker Genes'; end
if nargin < 3, Tm = []; end
if nargin < 2, sce = []; end

[missing, total] = pkg.i_missingmarkers(sce, Tm);
if isempty(missing), return; end

shown = missing;
if numel(shown) > 12
    shown = [shown(1:12); "..."];
end
gui.myWarndlg(parentfig, sprintf(['%d of %d marker gene(s) are not in ' ...
    'this dataset and will be ignored: %s. A cell type whose markers are ' ...
    'all missing can never be assigned, so check the gene symbols and the ' ...
    'species.'], numel(missing), total, strjoin(shown, ', ')), ...
    dlgtitle);
end
