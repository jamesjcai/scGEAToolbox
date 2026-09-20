function callback_CausalCCCNetworkView(src, ~)
%GUI.CALLBACK_CAUSALCCCNETWORKVIEW  Menu callback: pick sender and receiver
%cell groups, then plot a local causalCCC-style network for them.
%
%   Wires GUI.E_SCE2CAUSALCCCVIEW into the app's Network menu, next to the
%   other Cell-Cell Communication entries (One-/Two-Sample scTenifoldXct,
%   scTenifoldCko, Talklr).

[FigureHandle, sce] = gui.gui_getfigsce(src);

answer = gui.myQuestdlg(FigureHandle, sprintf(...
    ['This builds a LOCAL approximation of a causalCCC/MIIC network - ', ...
     'mutual-information (or correlation) edges within each side plus the ', ...
     'known ligand-receptor links, not the server''s causal discovery. ', ...
     'Use "Export causalCCC/MIIC Input Files..." for the real analysis on ', ...
     'https://miic.curie.fr/causalCCC.php. Continue?']));
if ~strcmp(answer, 'Yes'), return; end

% Cell type when it is there, cluster id when it is not, a warning when
% neither can supply two groups - see GUI.I_GETCELLGROUPS, which raises the
% dialogs and returns empty once the user has been told.
[labels, groupBy] = gui.i_getcellgroups(sce, FigureHandle);
if isempty(labels), return; end
cL = unique(labels(labels ~= ""), 'stable');

senders = i_selectonegroup(cL, 'Select the SENDER cell group:', FigureHandle);
if isempty(senders), return; end

remaining = cL(~ismember(cL, senders));
if isempty(remaining)
    gui.myWarndlg(FigureHandle, ...
        'No cell group is left to pick as receiver - not every group can be a sender.');
    return;
end
receivers = i_selectonegroup(remaining, 'Select the RECEIVER cell group:', FigureHandle);
if isempty(receivers), return; end

gui.e_sce2causalcccview(sce, senders, receivers, FigureHandle, 'GroupBy', groupBy);
end


function sel = i_selectonegroup(items, promptstr, parentfig)
% One item from a plain string list, single-selection only - the same
% lightweight listdlg/myListdlg pattern gui.i_select1class uses, instead of
% gui.i_selmultidialog's two-pane shuttle-box table, which is built for
% picking several items at once and is unnecessarily heavy (and misleading
% - it invites multi-selection) for choosing a single sender or receiver
% group here.
sel = strings(1, 0);
if gui.i_isuifig(parentfig)
    [indx, tf] = gui.myListdlg(parentfig, cellstr(items), promptstr, [], false);
else
    [indx, tf] = listdlg('PromptString', {promptstr}, ...
        'SelectionMode', 'single', 'ListString', cellstr(items), ...
        'ListSize', [220, 300], 'Name', ' ');
end
if tf == 1
    sel = items(indx);
end
end
