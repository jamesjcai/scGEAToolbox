function [s] = e_fgsearun(T, rmribo)
%E_FGSEARUN  Preranked gene set enrichment for a scTenifold gene ranking.
%
%   s = TEN.E_FGSEARUN(T) tests the ranked gene list in T.genelist against
%   the Enrichr libraries and returns the enriched sets. T is the table
%   scTenifoldNet and scTenifoldKnk produce, already ordered by regulatory
%   distance, so the order of T.genelist is the ranking.
%
%   Native MATLAB by way of SC_FGSEA. The R fgsea package is still reachable
%   as RUN.R_FGSEA if a like-for-like comparison is wanted.
%
%   INPUTS:
%     T      - table with a genelist column, most interesting gene first.
%     rmribo - drop ribosomal genes before ranking (default true).
%
%   See also SC_FGSEA, SC_GSETTEST, RUN.R_FGSEA.

if nargin < 2, rmribo = true; end

if istable(T)
    if ~ismember('genelist', T.Properties.VariableNames)
        error("ten:e_fgsearun:noGenelist", ...
            "T needs a genelist column; it has %s.", ...
            strjoin(T.Properties.VariableNames, ", "));
    end
    genelist = T.genelist;
else
    genelist = T;
end

s = sc_fgsea(genelist, [], RemoveRibosomal=rmribo);
end
