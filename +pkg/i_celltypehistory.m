function [names, labels] = i_celltypehistory(sce)
%I_CELLTYPEHISTORY Every cell type annotation held on SCE, newest first.
%
%   [names, labels] = pkg.i_celltypehistory(sce)
%
%   NAMES names each annotation: "Current" for the active c_cell_type_tx,
%   then each 'old_cell_type_N' cell attribute that PKG.I_STASHCELLTYPEHISTORY
%   wrote, highest N first so the most recently replaced annotation follows
%   the active one. LABELS is a cell array of string arrays, one per name,
%   each with sce.NumCells elements.
%
%   The active annotation is listed even when it is all "undetermined", so a
%   caller can show that an annotation was replaced by nothing rather than
%   silently omitting the current state.
%
%   See also PKG.I_STASHCELLTYPEHISTORY, GUI.CALLBACK_COMPARECELLTYPEANNOTATIONS.

names = string.empty(1, 0);
labels = {};
if nargin < 1 || isempty(sce) || sce.NumCells == 0, return; end

names = "Current";
labels = {string(sce.c_cell_type_tx(:))};

attrnames = string(sce.list_cell_attributes(1:2:end));
num = nan(1, numel(attrnames));
for k = 1:numel(attrnames)
    tok = regexp(attrnames(k), '^old_cell_type_(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        num(k) = str2double(tok{1});
    end
end

% Highest N is the most recently stashed, so descending order puts the
% annotations in the order they were replaced, newest first.
found = find(~isnan(num));
[~, ord] = sort(num(found), 'descend');
for k = found(ord)
    v = string(sce.list_cell_attributes{2*k});
    if numel(v) ~= sce.NumCells, continue; end
    names(end+1) = attrnames(k);   %#ok<AGROW>
    labels{end+1} = v(:);          %#ok<AGROW>
end
end
