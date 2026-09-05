function tf = i_confirmoverwritecelltype(FigureHandle, sce)
%I_CONFIRMOVERWRITECELLTYPE Ask before an annotation method replaces cell types.
% Returns true if it is fine to proceed: either there is nothing worth
% keeping, or the user confirmed. Every caller - SC_ANNOTATECELLS, the
% reference-model callbacks (GUI.CALLBACK_RUNSCIMILARITY,
% GUI.CALLBACK_RUNPANHUMANPY), marker-gene annotation
% (GUI.CALLBACK_DETERMINECELLTYPECLUSTERS),
% GUI.CALLBACK_ASSIGNCELLTYPEFROMATTRIB and GUI.SC_CELLATTRIBEDITOR -
% stashes the current labels into a new 'old_cell_type_N' cell attribute
% (see PKG.I_STASHCELLTYPEHISTORY) before overwriting, so this is about the
% active c_cell_type_tx changing, not about data loss.
%
% SINGLECELLEXPERIMENT.ASSIGNCELLTYPE stashes too, but has no figure to ask
% in; its keepold flag stands in for this dialog.

tf = true;
if ~pkg.i_hascelltypelabels(sce), return; end

msg = ['This replaces the current, active cell type labels. They will be ', ...
    'kept as a new ''old_cell_type_N'' cell attribute, which Annotate > ', ...
    'Compare Cell Type Annotations shows beside the new labels and ', ...
    'Annotate > Assign Cell Type from Cell Attribute Table switches back ', ...
    'to. Continue?'];
tf = strcmp(gui.myQuestdlg(FigureHandle, msg, 'Annotate Cell Types'), 'Yes');
end
