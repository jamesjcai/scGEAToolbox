function tf = i_confirmoverwritecelltype(FigureHandle, sce)
%I_CONFIRMOVERWRITECELLTYPE Ask before an annotation method replaces cell types.
% Returns true if it is fine to proceed: either there is nothing worth
% keeping, or the user confirmed. SC_ANNOTATECELLS and the reference-model
% callbacks (GUI.CALLBACK_RUNSCIMILARITY, GUI.CALLBACK_RUNPANHUMANPY) all
% stash the current labels into a new 'old_cell_type_N' cell attribute (see
% PKG.I_STASHCELLTYPEHISTORY) before overwriting, so this is about the
% active c_cell_type_tx changing, not about data loss.

tf = true;
if ~pkg.i_hascelltypelabels(sce), return; end

msg = ['This replaces the current, active cell type labels. The current ', ...
    'labels will be kept as a new ''old_cell_type_N'' cell attribute ', ...
    '(View/Export Cell Attribute Table to see them, or Assign Cell Type ', ...
    'from Cell Attribute Table to restore them). Continue?'];
tf = strcmp(gui.myQuestdlg(FigureHandle, msg, 'Annotate Cell Types'), 'Yes');
end
