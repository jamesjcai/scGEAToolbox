function [needupdatesce] = callback_AnnotateByReference(src, ~, source)
%CALLBACK_ANNOTATEBYREFERENCE Label clusters against a reference the user supplies.
%
%   needupdatesce = gui.callback_AnnotateByReference(app)
%   needupdatesce = gui.callback_AnnotateByReference(app, [], "allen")
%
% Annotate > Atlas-Based Cell Type Annotation > Custom Reference, and with
% SOURCE "allen" the Allen Brain Cell Atlas item beside it. The
% reference is either an annotated SCE, which is summarised on the spot by
% SC_BUILDCELLTYPEREF, or a reference struct already made by it (or by hand
% from a published atlas's tables), from a .mat file or the workspace.
% Labelling goes through SC_ANNOTATECELLS(Method="reference"), so the old
% labels are stashed and the reference is recorded in sce.metadata exactly
% as on the command line.
%
% Returns true when c_cell_type_tx was replaced and the app must refresh.
%
% see also: sc_celltypeannoref, sc_buildcelltyperef, sc_annotatecells

if nargin < 3, source = ""; end
needupdatesce = false;
[FigureHandle, sce] = gui.gui_getfigsce(src);

% The method labels clusters, so it has nothing to work on without them.
% SingleCellExperiment seeds c_cluster_id to all-ones.
if isempty(sce.c_cluster_id) || isscalar(unique(string(sce.c_cluster_id)))
    gui.myWarndlg(FigureHandle, ['Reference annotation labels clusters, ' ...
        'and the cells have not been clustered. Cluster them first ' ...
        '(Cluster > Cluster Cells).']);
    return;
end

ref = gui.i_pickcelltyperef(FigureHandle, source, sce.g);
if isempty(ref), return; end

if ~gui.i_confirmoverwritecelltype(FigureHandle, sce), return; end

fw = gui.myWaitbar(FigureHandle);
try
    [~, ~, stashname, Tref] = sc_annotatecells(sce, Method="reference", ...
        Reference=ref, Verbose=false);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);
needupdatesce = true;

numDisagree = sum(~Tref.Agree);
msg = sprintf('%d clusters labelled against %s.', height(Tref), ref.Source);
if numDisagree > 0
    msg = sprintf(['%s The correlation and marker tests disagree on %d ' ...
        'of them; see the rows with Agree = false in the table.'], ...
        msg, numDisagree);
end
gui.myHelpdlg(FigureHandle, msg + gui.i_stashnotice(stashname));
gui.TableViewerApp(Tref, FigureHandle, "ReferenceAnnotation");
end
