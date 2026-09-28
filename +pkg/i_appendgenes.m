function [X, g] = i_appendgenes(X, g, extrag, Xall, gall)
%I_APPENDGENES Add genes to a gene subset, taking counts from the full matrix.
%
%   [X, g] = PKG.I_APPENDGENES(X, g, extrag, Xall, gall) appends the
%   entries of EXTRAG that are not already in G to the subset (X, g),
%   reading their counts out of the full matrix (XALL, GALL). Genes already
%   present in G are left where they are, so an HVG ranking keeps its order
%   and only gains rows at the bottom.
%
%   Used to force a whitelist - PanglaoDB markers, say - into a gene set
%   that was selected on variability alone. Canonical markers of rare cell
%   types are often expressed in too few cells to rank as highly variable,
%   so they drop out exactly when they are needed.
%
%   EXTRAG must be a subset of GALL; anything else is a caller bug and
%   raises rather than silently selecting the wrong rows.
%
%   See also PKG.I_GETMARKERWHITELIST.

g = string(g(:));
gall = string(gall(:));
extrag = unique(string(extrag(:)));

missing = extrag(~ismember(extrag, g));
if isempty(missing), return; end

[tf, loc] = ismember(missing, gall);
if ~all(tf)
    error('pkg:i_appendgenes:geneNotFound', ...
        ['%d of the genes to append are not in GALL, the first being ' ...
        '"%s". Intersect the list with the gene list first.'], ...
        sum(~tf), missing(find(~tf, 1)));
end

X = [X; Xall(loc, :)];
g = [g; missing];
end
