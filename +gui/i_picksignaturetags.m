function [rows, taglabel] = i_picksignaturetags(T, parentfig)
%I_PICKSIGNATURETAGS Narrow the predefined score table to one or more collections.
%
%   [rows, taglabel] = gui.i_picksignaturetags(T, parentfig)
%
%   Input:
%     T         - the table returned by pkg.e_cellscores(X, g, 0). Its
%                 SignatureTag column names the collection each signature
%                 came from: CancerSEA, ISGs, Glycobiology, Cell Cycle
%                 [MSigDB] and a dozen more.
%     parentfig - figure the dialog should be parented to. Optional.
%
%   Outputs:
%     rows     - row indices into T, in table order; empty if the user
%                cancelled
%     taglabel - the chosen collection(s), for the next dialog's prompt, or
%                "" when the whole table was chosen
%
%   Without this step a score dialog is one flat list of every signature in
%   the table, so a user who knows which collection they want still has to
%   scroll past all the others. Picking "All signatures" keeps that flat
%   list, so the extra step costs one click for anyone browsing everything.
%
%   Shared by gui.callback_CompareCellScoreBtwCls and
%   gui.callback_GetCellSignatureMatrix so both present the same dialog.

if nargin < 2, parentfig = []; end

rows = [];
taglabel = "";

tags = string(T.SignatureTag);
% Rows the spreadsheet leaves untagged still have to be reachable.
tags(ismissing(tags) | strlength(strtrim(tags)) == 0) = "(untagged)";

[uTags, ~, tagid] = unique(tags);          % unique() sorts, which is the
counts = accumarray(tagid, 1);             % order we want in the dialog

allLabel = "All signatures (" + string(height(T)) + ")";
listitems = [allLabel; uTags + " (" + string(counts) + ")"];

if gui.i_isuifig(parentfig)
    [indx, tf] = gui.myListdlg(parentfig, listitems, 'Select Collection:');
else
    [indx, tf] = listdlg('PromptString', 'Select Collection', ...
        'SelectionMode', 'multiple', 'ListString', ...
        cellstr(listitems), 'ListSize', [300, 300]);
end
if tf ~= 1 || isempty(indx), return; end

% Index 1 is the "All signatures" row; the rest are offset by it.
if any(indx == 1)
    rows = (1:height(T))';
    return;
end

chosen = uTags(indx - 1);
rows = find(ismember(tags, chosen));
% Tag names run long (e.g. "CancerSEA [PMID:30329142]"), so naming more
% than a couple of them would push the next dialog's prompt off the edge.
if numel(chosen) <= 2
    taglabel = strjoin(chosen(:)', ', ');
else
    taglabel = string(numel(chosen)) + " collections";
end
end
