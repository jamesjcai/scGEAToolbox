function name = i_stashcelltypehistory(sce)
%I_STASHCELLTYPEHISTORY Save the current cell type labels before overwriting.
% Adds a new 'old_cell_type_N' cell attribute (SCE is a handle, so this
% mutates it in place) rather than replacing a single 'old_cell_type' slot,
% so every past annotation survives instead of just the most recent one.
% Returns the attribute name used, or "" if there were no labels worth
% keeping (see I_HASCELLTYPELABELS).

name = "";
if ~pkg.i_hascelltypelabels(sce), return; end

existing = string(sce.list_cell_attributes(1:2:end));
n = 0;
for k = 1:numel(existing)
    tok = regexp(existing(k), '^old_cell_type_(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        n = max(n, str2double(tok{1}));
    end
end
name = sprintf('old_cell_type_%d', n + 1);
sce.setCellAttribute(name, string(sce.c_cell_type_tx));
end
