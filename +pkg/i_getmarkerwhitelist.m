function [whitelist] = i_getmarkerwhitelist(g, X)
%I_GETMARKERWHITELIST PanglaoDB marker genes present in a gene list.
%
%   whitelist = PKG.I_GETMARKERWHITELIST(g) returns the entries of G that
%   are PanglaoDB cell-type markers, in G's own spelling. Matching is
%   case-insensitive, so the uppercase reference list ("CD3E") also picks
%   up mouse symbols ("Cd3e").
%
%   whitelist = PKG.I_GETMARKERWHITELIST(g, X) additionally drops markers
%   that carry no counts in any cell. A whitelist exists to add genes back
%   to an HVG ranking that SC_SPLINEFIT built after removing all-zero
%   genes; putting an all-zero gene back only widens the matrix.
%
%   These are the ~4700 markers of PRIMARY cell types. For subdividing one
%   type into its subtypes, the subtype markers of that type are the
%   relevant list and this one is mostly noise - see SC_CSUBTYPEANNO, which
%   matches its own table with PKG.I_MATCHGENENAMES instead.
%
%   See also PKG.I_GET_PANGLAODBMARKERS, PKG.I_MATCHGENENAMES.

if nargin < 2, X = []; end

whitelist = strings(0, 1);
markerg = pkg.i_get_panglaodbmarkers;

% I_GET_PANGLAODBMARKERS swallows a failed load and hands back "".
if isempty(markerg) || (isscalar(markerg) && markerg == "")
    return;
end

whitelist = pkg.i_matchgenenames(g, markerg, X);
end
