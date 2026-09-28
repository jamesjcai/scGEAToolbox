function [es, pk] = e_gseascore(w, q, n, m)
%E_GSEASCORE  Weighted Kolmogorov-Smirnov enrichment score for one gene set.
%
%   es = PKG.E_GSEASCORE(w, q, n, m) returns the GSEA running-sum
%   enrichment score of one gene set against one ranked gene list.
%
%   [es, pk] = PKG.E_GSEASCORE(___) also returns the position in the ranked
%   list at which the running sum reaches its extremum, which is the end of
%   the leading edge.
%
%   THE POINT OF THIS FORM. The running sum only changes slope at a hit, so
%   its extremum can only land on one -- at a hit for a positive score, and
%   immediately before one for a negative score. Evaluating it at the M hit
%   positions rather than at all N ranks makes each score O(M log M) instead
%   of O(N). That is not a micro-optimisation: every permutation null needs
%   one score per set per draw, so with thousands of sets and thousands of
%   draws the O(N) form is the difference between a run of minutes and a run
%   that is abandoned.
%
%   INPUTS:
%     w - N-by-1 weights in ranked order, abs(stat).^p for the usual
%         weighted score, or all ones for the classic unweighted one.
%     q - M-by-1 ASCENDING positions of the set's genes within that ranking.
%         Ascending is not checked; passing unsorted positions silently
%         returns a wrong score. Taking a column of a sparse membership
%         matrix with FIND gives them in the right order for free.
%     n - length of the ranked list.
%     m - number of set genes in it, numel(q).
%
%   OUTPUTS:
%     es - enrichment score in [-1, 1], the signed maximum deviation of the
%          running sum from zero. Positive means the set sits toward the
%          top of the ranking.
%     pk - rank at which that deviation occurs, 1-based. 0 when ES is 0,
%          i.e. when there is no leading edge to report.
%
%   Callers turn PK into a leading edge by taking the set genes at ranks
%   <= PK when ES > 0, and at ranks >= PK when ES < 0.
%
% REF: Subramanian et al. (2005) PNAS 102:15545.
%
% See also RUN.ML_GSEA, SC_GSETTEST, SC_FGSEA.

nr = sum(w(q));
if nr <= 0 || m >= n || m < 1
    es = 0;
    pk = 0;
    return;
end
cumhit = cumsum(w(q))/nr;
misses = (q - (1:m)')/(n - m);

% Both candidates are padded with a leading zero, so that a set whose
% running sum never leaves one side of zero reports that side as 0 rather
% than as its least extreme excursion. The padding shifts the index by one.
[espos, ip] = max([0; cumhit - misses]);            % maxima land on a hit
[esneg, in] = min([0; [0; cumhit(1:m-1)] - misses]); % minima just before one

if espos > -esneg
    es = espos;
    if ip > 1
        pk = q(ip - 1);
    else
        pk = 0;
    end
else
    es = esneg;
    if in > 1
        pk = q(in - 1) - 1;
    else
        pk = 0;
    end
end
if es == 0
    pk = 0;
end

end
