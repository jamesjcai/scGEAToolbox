function [done] = callback_HarmonyAny(src, ~)
%CALLBACK_HARMONYANY Batch integration with Harmony, backend chosen once.
%
%   Harmony used to appear three times in the menus -- under Edit as the
%   native MATLAB run, and again under R Tools and Python Tools -- with the
%   same label on all three. A user who found one had no way to know the
%   other two existed, or which to prefer. This is the single entry point:
%   it asks which backend to use, offering R and Python only when they are
%   configured, and hands off to the callback that runs it.
%
%   The native implementation is the default because RUN.ML_HARMONY needs
%   no installation, and because it is the only one of the three that can
%   correct the principal components and re-embed rather than correcting
%   the two dimensions already on screen.
%
%   SRC is the app, so the orchestration the app menu used to do around
%   these calls -- colouring by batch, refreshing, offering to store the
%   corrected embedding -- happens here too.
%
%   See also GUI.CALLBACK_HARMONY, GUI.CALLBACK_HARMONYR,
%   GUI.CALLBACK_HARMONYPY, RUN.ML_HARMONY.

done = false;
[FigureHandle, sce] = gui.gui_getfigsce(src);

if ~gui.gui_showrefinfo('Harmony [PMID:31740819]', FigureHandle), return; end

if numel(unique(sce.c_batch_id)) < 2
    gui.myWarndlg(FigureHandle, ...
        'No batch effect (all cells have the same SCE.C_BATCH_ID)');
    return;
end

% Only offer a backend that can actually run. Both alternatives fall back
% to the native implementation when their runtime is missing, so an
% unconfigured entry would be a choice that silently does something else.
options = {'MATLAB (native)'};
hasR = ispref('scgeatoolbox', 'rexecutablepath') && ...
    ~isempty(getpref('scgeatoolbox', 'rexecutablepath', []));
if hasR
    options{end+1} = 'R (harmony package)';
end
if pkg.i_checkpython
    options{end+1} = 'Python (harmonypy)';
end

backend = options{1};
if numel(options) > 1
    backend = gui.myQuestdlg(FigureHandle, 'Choose Harmony backend:', ...
        'Harmony', options, options{1});
    if isempty(backend), return; end
end

switch backend
    case 'R (harmony package)'
        fun = @gui.callback_HarmonyR;
    case 'Python (harmonypy)'
        fun = @gui.callback_Harmonypy;
    otherwise
        fun = @gui.callback_Harmony;
end

% Harmony corrects the embedding, so the plot is only readable afterwards
% if the cells are coloured by the batch it corrected for.
if isa(src, 'matlab.apps.AppBase')
    c1 = findgroups(string(src.sce.c));
    c2 = findgroups(string(src.sce.c_batch_id));
    if ~isequal(c1, c2)
        answer = gui.myQuestdlg(FigureHandle, ...
            'Color cells by batch id (SCE.C_BATCH_ID)?', '');
        switch answer
            case 'Yes'
                [src.c, src.cL] = findgroups(string(src.sce.c_batch_id));
                src.sce.c = src.c;
                src.in_RefreshAll(true, false);
            case 'No'
            otherwise
                return;
        end
    end
end

pcsBefore = i_harmonypcs(sce);
if ~fun(src), return; end
done = true;

if ~isa(src, 'matlab.apps.AppBase'), return; end
[src.c, src.cL] = findgroups(string(src.sce.c));
src.in_RefreshAll(true, false);

% Only the native backend's "Correct PCs" path leaves corrected components,
% and the existing clusters were computed before them. Offer to redo the
% clustering on them now, while it is obvious why.
pcsAfter = i_harmonypcs(src.sce);
% New corrected PCs mean the native "Correct PCs" path ran, and that path
% always redraws the cells as a UMAP (GUI.CALLBACK_HARMONY).
reembedded = ~isempty(pcsAfter) && ~isequal(pcsBefore, pcsAfter);
if reembedded
    if strcmp('Yes', gui.myQuestdlg(FigureHandle, ...
            ['Re-cluster cells on the batch-corrected principal ' ...
            'components (Louvain, resolution 0.8)?'], ''))
        fw = gui.myWaitbar(FigureHandle);
        try
            src.sce.clustercells([], 'louvainpc', true);
        catch ME
            gui.myWaitbar(FigureHandle, fw, true);
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end
        gui.myWaitbar(FigureHandle, fw);
        [src.c, src.cL] = findgroups(string(src.sce.c_cluster_id));
        src.sce.c = src.c;
        src.in_RefreshAll(true, false);
    end
end

dim = size(src.sce.s, 2);
if reembedded
    q = sprintf('Save the corrected UMAP as the stored UMAP %dD embedding?', dim);
else
    q = 'Update Saved Embedding?';
end
if ~strcmp('Yes', gui.myQuestdlg(FigureHandle, q, '')), return; end
if reembedded
    % A UMAP goes in the UMAP slot only. The picker below offered tSNE and
    % PHATE too, so a UMAP could be stored under tsne2d.
    methodtag = sprintf('umap%dd', dim);
else
    [methodtag] = gui.i_pickembedmethod(FigureHandle, false, dim);
    if isempty(methodtag), return; end
    if iscell(methodtag), methodtag = cell2mat(methodtag); end
end
if ismember(methodtag, fieldnames(src.sce.struct_cell_embeddings))
    src.sce.struct_cell_embeddings.(methodtag) = src.sce.s;
    gui.myHelpdlg(FigureHandle, ...
        sprintf('%s Embedding is updated.', methodtag));
end

end

function pcs = i_harmonypcs(sce)
pcs = [];
if isfield(sce.struct_cell_reductions, 'harmony')
    pcs = sce.struct_cell_reductions.harmony;
end
end
