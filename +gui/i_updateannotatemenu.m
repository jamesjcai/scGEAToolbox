function i_updateannotatemenu(app)
%I_UPDATEANNOTATEMENU Set the Annotate menu's item states as it opens.
%
%   gui.i_updateannotatemenu(app)
%
% Called from the Annotate menu's own MenuSelectedFcn, which MATLAB runs when
% the menu is expanded. Only Enable is touched: adding or removing items from a
% menu callback can leave the menu blank (see the uimenu documentation).
%
% Right now this greys out the two subtype annotation items - the bundled one
% unless the data holds a cell type the bundled subtype marker table has
% markers for, the customized one unless the cells are annotated at all - and
% says why in the tooltip. Anything else that should depend on what is in the
% data belongs here too.
%
% see also: pkg.i_subtypecandidates, gui.callback_SubtypeAnnotation

if nargin < 1 || isempty(app), return; end

sce = [];
if isprop(app, 'sce'), sce = app.sce; end

hascelltypes = ~isempty(sce) && isa(sce, 'SingleCellExperiment') && ...
    pkg.i_hascelltypelabels(sce);
[ctypelist, ~, primarytypes] = pkg.i_subtypecandidates(sce);

in_setbundled(app, hascelltypes, ctypelist, primarytypes);
in_setcustom(app, hascelltypes);
end

function in_setbundled(app, hascelltypes, ctypelist, primarytypes)
if ~isprop(app, 'AnnotateCellSubtypesMenu'), return; end
m = app.AnnotateCellSubtypesMenu;
if ~pkg.i_isvalid(m), return; end

if ~isempty(ctypelist)
    m.Enable = 'on';
    m.Tooltip = {sprintf('Re-cluster and label the subtypes of: %s', ...
        strjoin(ctypelist, ', '))};
    return;
end

% Two different reasons to be greyed out, and the fix differs: annotate the
% cells, or use the customized marker item just below.
m.Enable = 'off';
if ~hascelltypes
    m.Tooltip = {sprintf(['Cells are not annotated yet. Annotate cell ' ...
        'types first. Subtype markers are available for: %s.'], ...
        strjoin(primarytypes, ', '))};
else
    m.Tooltip = {sprintf(['None of the cell types in this dataset has ' ...
        'bundled subtype markers (available for: %s). Use customized ' ...
        'marker genes instead.'], strjoin(primarytypes, ', '))};
end
end

function in_setcustom(app, hascelltypes)
% The customized item asks nothing of the marker table, so annotated cells are
% the whole requirement.

if ~isprop(app, 'AnnotateCellSubtypesUsingCustomizedMarkersMenu'), return; end
m = app.AnnotateCellSubtypesUsingCustomizedMarkersMenu;
if ~pkg.i_isvalid(m), return; end

if hascelltypes
    m.Enable = 'on';
    m.Tooltip = {['Re-cluster any annotated cell type and label its ' ...
        'subtypes against a marker list you supply']};
else
    m.Enable = 'off';
    m.Tooltip = {['Cells are not annotated yet. Annotate cell types ' ...
        'first, then come back for their subtypes.']};
end
end
