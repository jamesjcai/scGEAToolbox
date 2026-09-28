function [matched] = i_matchgenenames(g, wanted, X)
%I_MATCHGENENAMES Genes of a list named by another list, in the list's spelling.
%
%   matched = PKG.I_MATCHGENENAMES(g, wanted) returns the entries of G that
%   appear in WANTED, in G's own spelling. Matching is case-insensitive, so
%   an uppercase reference list ("CD3E") also picks up mouse symbols
%   ("Cd3e").
%
%   Returning G's spelling rather than WANTED's is the point of the
%   function: the result is handed to code that indexes back into G, and an
%   uppercased name matches nothing in a mouse dataset. Marker tables in
%   this toolbox are upper-cased by convention - PKG.E_MARKERWEIGHT does it
%   - so every caller mixing a marker list with a gene list needs this.
%
%   matched = PKG.I_MATCHGENENAMES(g, wanted, X) additionally drops genes
%   that carry no counts in any cell, where X has one row per gene in G.
%
%   See also PKG.I_GETMARKERWHITELIST, PKG.I_APPENDGENES.

if nargin < 3, X = []; end

matched = strings(0, 1);
if isempty(wanted), return; end

g = string(g(:));
wanted = unique(upper(string(wanted(:))));
wanted(wanted == "") = [];
if isempty(wanted), return; end

idx = ismember(upper(g), wanted);

if ~isempty(X)
    assert(size(X, 1) == numel(g), ...
        'PKG.I_MATCHGENENAMES: X must have one row per gene in G.');
    idx(idx) = full(sum(X(idx, :) ~= 0, 2)) > 0;
end

matched = unique(g(idx));
end
