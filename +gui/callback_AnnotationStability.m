function [needupdatesce, speciestag] = callback_AnnotationStability(src, ~)
%CALLBACK_ANNOTATIONSTABILITY Pick a clustering resolution by how its labels behave.
%
%   [needupdatesce, speciestag] = gui.callback_AnnotationStability(app)
%
% Cluster > Choose Clustering Resolution by Annotation. Runs
% SC_ANNOTATIONSTABILITY through GUI.I_ANNOTATIONSWEEP, which asks for the
% labelling method (database markers, a custom marker list or a reference),
% runs the sweep once per dataset and offers the saved result after that,
% and shows the per-resolution table with its Methods Text. It then offers
% to apply the recommendation: recluster at that resolution, label the
% clusters by the same method through SC_ANNOTATECELLS (which stashes the
% old labels and records provenance), and keep each cell's stability
% across resolutions as the 'annotation_stability' cell attribute.
%
% Cluster Cells > Louvain on PCs > Choose by Annotation
% (GUI.CALLBACK_RECLUSTERCELLS) uses the same sweep but sets only the
% clusters.
%
% Returns true when clusters and labels were replaced and the app must
% refresh, and the species chosen for the marker method ([] otherwise).
%
% see also: sc_annotationstability, gui.i_annotationsweep, gui.i_pickcelltyperef

needupdatesce = false;
[FigureHandle, sce] = gui.gui_getfigsce(src);

[S, args, speciestag] = gui.i_annotationsweep(src, FigureHandle, sce);
if isempty(S), return; end

best = S.T(S.T.Recommended, :);
msg = sprintf(['Recommended resolution %.3g: %d clusters, %d cell types. ' ...
    'The table shows every resolution tried.\n\nApply reclusters the cells ' ...
    'at this resolution and annotates the clusters. Not Now leaves the ' ...
    'clusters and cell types as they are; the sweep is saved with the ' ...
    'data, so running this again offers it without recomputing.'], ...
    S.Resolution, best.NumClusters, best.NumTypes);
if ~strcmp(gui.myQuestdlg(FigureHandle, msg, 'Apply Resolution', ...
        {'Apply', 'Not Now'}, 'Apply'), 'Apply')
    return;
end
if S.Method == "reference" && isempty(args)
    % A saved reference sweep: the reference itself is not kept with the
    % data, so ask for it again to label the clusters.
    gui.myHelpdlg(FigureHandle, sprintf(['Select the reference the sweep ' ...
        'used again (%s) to label the clusters.'], S.MethodLabel));
    ref = gui.i_pickcelltyperef(FigureHandle, "", sce.g);
    if isempty(ref), return; end
    args = {'Reference', ref};
end
if ~gui.i_confirmoverwritecelltype(FigureHandle, sce), return; end

fw = gui.myWaitbar(FigureHandle);
try
    sce.c_cluster_id = S.ClusterId;
    sce.struct_cell_clusterings.louvainpc = S.ClusterId;
    [~, ~, stashname] = sc_annotatecells(sce, 'Method', S.Method, args{:}, ...
        'Verbose', false);
    sce.setCellAttribute('annotation_stability', S.Stability);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);
needupdatesce = true;

gui.myHelpdlg(FigureHandle, sprintf(['Reclustered at resolution %.3g ' ...
    'and annotated. Each cell''s share of resolutions that agreed on its ' ...
    'type is kept as the ''annotation_stability'' cell attribute; low ' ...
    'values mark cells whose type depends on the resolution.'], ...
    S.Resolution) + gui.i_stashnotice(stashname));
end
