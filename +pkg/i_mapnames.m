function out = i_mapnames(names, vars, cols)
%PKG.I_MAPNAMES  Translate gene symbols into the dataset/node column names
%they were given after PKG.I_CAUSALCCCASSEMBLE's disambiguation (a symbol
%picked as an intracellular gene on both sides is suffixed
%"_senders"/"_receivers"). Not in a "private" folder for the same reason as
%PKG.I_CAUSALCCCASSEMBLE: both SC_CAUSALCCCNET (top-level) and
%PKG.SC_SCE2CAUSALCCC call it.
%
%   out = PKG.I_MAPNAMES(names, vars, cols)
%
%   names - symbols to translate (e.g. a pair's ligand/receptor column)
%   vars  - the pre-disambiguation variable names (sndVars or rcvVars)
%   cols  - the matching post-disambiguation column names (sndCols/rcvCols)
[found, idx] = ismember(upper(names(:)), upper(vars(:)));
assert(all(found), ...
    'pkg.sc_sce2causalccc: %s missing from the exported columns.', ...
    strjoin(names(~found), ', '));
out = cols(idx);
out = out(:);
end
