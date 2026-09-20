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

if ~fun(src), return; end
done = true;

if ~isa(src, 'matlab.apps.AppBase'), return; end
[src.c, src.cL] = findgroups(string(src.sce.c));
src.in_RefreshAll(true, false);

if ~strcmp('Yes', gui.myQuestdlg(FigureHandle, ...
        'Update Saved Embedding?', '')), return; end
dim = size(src.sce.s, 2);
[methodtag] = gui.i_pickembedmethod(FigureHandle, false, dim);
if isempty(methodtag), return; end
if iscell(methodtag), methodtag = cell2mat(methodtag); end
if ismember(methodtag, fieldnames(src.sce.struct_cell_embeddings))
    src.sce.struct_cell_embeddings.(methodtag) = src.sce.s;
    gui.myHelpdlg(FigureHandle, ...
        sprintf('%s Embedding is updated.', methodtag));
end

end
