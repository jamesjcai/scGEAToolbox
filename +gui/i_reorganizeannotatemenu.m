function i_reorganizeannotatemenu(app)
%I_REORGANIZEANNOTATEMENU Move "choose method" next to the marker-based items.
% scgeatoolApp.mlapp is a binary App Designer file, so the Annotate menu is
% reorganized at launch time rather than by editing the .mlapp. Moves
% AnnotateCellTypesChooseMethodMenu (baked in as the last item, "Annotate Cell
% Types (choose method)...") up to just after "Annotate Cell Types Using
% Customized Marker Genes...", renamed "Annotate Cell Types (Other
% Methods)...", since GUI.CALLBACK_ANNOTATECELLS no longer offers a
% marker-gene-based option (that stays covered by the first two menu items).
% Idempotent: if the menu no longer exists at the expected place (e.g. later
% baked into the app in this new position), this does nothing.

if nargin < 1 || isempty(app), return; end
if ~isprop(app, 'AnnotateMenu') || ~pkg.i_isvalid(app.AnnotateMenu), return; end
requiredprops = ["AnnotateCellTypesChooseMethodMenu", ...
    "AnnotateCellTypesUsingPanglaoDBMenu", ...
    "AnnotateCellTypesUsingCustomizedMarkersMenu", ...
    "AssignCellTypefromCellAttributeTable", ...
    "ImportCellAnnotationfromSCEinWorkspaceMenu", ...
    "ImportCellAnnotationfromSCEDataFileMenu", ...
    "AnnotateCellTypesforBrushedCellsMenu", ...
    "FindMarkerGenesforBrushedCellsMenu", ...
    "MakeMarkerGeneHeatmapMenu", ...
    "EstimateCellCyclePhaseMenu", ...
    "EstimateDifferentiationPotencyPMID33244588Menu", ...
    "EstimateStemnessPMID29625051Menu", ...
    "EstimateDissociationGeneRatioPMID34020534Menu", ...
    "SingleClickSolutionAnnotationMenu"];
if ~all(arrayfun(@(p) isprop(app, p), requiredprops)), return; end

choose = app.AnnotateCellTypesChooseMethodMenu;
choose.Text = 'Annotate Cell Types (Other Methods)...';
choose.Tooltip = {'LLM, SCimilarity, panhumanpy, or a consensus of them'};
choose.Separator = 'on';
app.AssignCellTypefromCellAttributeTable.Separator = 'on';

% Children is stored in reverse of display order, so build the desired
% top-to-bottom order and assign its flip.
displayorder = [ ...
    app.AnnotateCellTypesUsingPanglaoDBMenu
    app.AnnotateCellTypesUsingCustomizedMarkersMenu
    choose
    app.AssignCellTypefromCellAttributeTable
    app.ImportCellAnnotationfromSCEinWorkspaceMenu
    app.ImportCellAnnotationfromSCEDataFileMenu
    app.AnnotateCellTypesforBrushedCellsMenu
    app.FindMarkerGenesforBrushedCellsMenu
    app.MakeMarkerGeneHeatmapMenu
    app.EstimateCellCyclePhaseMenu
    app.EstimateDifferentiationPotencyPMID33244588Menu
    app.EstimateStemnessPMID29625051Menu
    app.EstimateDissociationGeneRatioPMID34020534Menu
    app.SingleClickSolutionAnnotationMenu];
app.AnnotateMenu.Children = flipud(displayorder);
end
