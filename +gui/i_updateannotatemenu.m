function i_updateannotatemenu(app)
%I_UPDATEANNOTATEMENU Set the Annotate menu's item states as it opens.
%
%   gui.i_updateannotatemenu(app)
%
% Called from the Annotate menu's own MenuSelectedFcn, which MATLAB runs when
% the menu is expanded. Only Enable is touched: adding or removing items from a
% menu callback can leave the menu blank (see the uimenu documentation).
%
% Right now this greys out "Annotate Cell Subtypes..." unless the data holds a
% cell type the bundled subtype marker table has markers for, and says why in
% the tooltip. Anything else that should depend on what is in the data belongs
% here too.
%
% see also: pkg.i_subtypecandidates, gui.callback_SubtypeAnnotation

if nargin < 1 || isempty(app), return; end
if ~isprop(app, 'AnnotateCellSubtypesMenu'), return; end

m = app.AnnotateCellSubtypesMenu;
if ~pkg.i_isvalid(m), return; end

sce = [];
if isprop(app, 'sce'), sce = app.sce; end

[ctypelist, ~, primarytypes] = pkg.i_subtypecandidates(sce);

if isempty(ctypelist)
    % Two different reasons to be greyed out, and the fix differs: annotate the
    % cells, or accept that this dataset holds no type the table covers.
    m.Enable = 'off';
    if isempty(sce) || ~pkg.i_hascelltypelabels(sce)
        m.Tooltip = {sprintf(['Cells are not annotated yet. Annotate cell ' ...
            'types first. Subtype markers are available for: %s.'], ...
            strjoin(primarytypes, ', '))};
    else
        m.Tooltip = {sprintf(['None of the cell types in this dataset has ' ...
            'subtype markers. Available for: %s.'], ...
            strjoin(primarytypes, ', '))};
    end
    return;
end

m.Enable = 'on';
m.Tooltip = {sprintf('Re-cluster and label the subtypes of: %s', ...
    strjoin(ctypelist, ', '))};
end
