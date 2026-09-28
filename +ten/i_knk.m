function [T] = i_knk(A0, idx, genelist, dosort, lambdav, reference)
% Virtual knockout of gene IDX in the network A0
%
% T = ten.i_knk(A0, idx, genelist)
% T = ten.i_knk(A0, idx, genelist, dosort, lambdav, reference)
%
% A0 is the network as TEN.SCTENIFOLDKNK returns it, in NET.PCRNET's
% orientation: row k holds the coefficients of the regression of gene k, so
% A0(k, j) is the effect of gene j ON gene k.
%
% THE KNOCKOUT REMOVES THE GENE'S OUTGOING EDGES (fixed 2026-09-25). A0 is
% transposed first, so that row k lists the genes k regulates, and that row
% is zeroed - as scTenifoldKnk in R does (WT <- t(WT); KO[gKO, ] <- 0). This
% used to zero row IDX of A0 itself, which removed the gene's regulators
% rather than its targets; and because TEN.I_MA symmetrizes the joint
% adjacency, the half-weight edges left at the gene were exactly the outgoing
% ones a knockout should remove. The link-count guard below now counts
% outgoing edges too.
%
% In practice the fix moves results little. The CP-denoised networks
% TEN.SCTENIFOLDKNK builds are nearly symmetric (||A - A'||/||A|| = 0.04-0.05
% on three 1000-gene populations), so the two perturbations nearly coincide:
% for the same knockout the two drdist rankings correlate at Spearman 0.998
% (median over 40 knockouts; the lowest 0.61), and every summary in
% pipelines/tenifold_laplacian/README.md ("Knockout direction") moved by at
% most a few percent. Results from before the fix can differ for individual
% knockouts, most for genes whose incoming and outgoing edges differ most.
%
% REFERENCE (default false) runs the rest as the R package does, for use by
% the reference mode of TEN.SCTENIFOLDKNK: strictDirection keeps ties,
% there is no 10-link minimum (only a warning when the gene has no outgoing
% edge), the alignment uses the unnormalized Laplacian on the unsymmetrized
% adjacency, and TEN.I_DR runs in its reference mode.
%
% see also: TEN.SCTENIFOLDKNK, TEN.KNKNULLCONTROL, TEN.I_MA, TEN.I_DR

if nargin < 6, reference = false; end
if nargin < 5, lambdav = 0; end
if nargin < 4, dosort = true; end
if ischar(idx) || isstring(idx)
    [~, idx] = ismember(idx, genelist);
end

import ten.*

if lambdav ~= 0
    if reference
        % strictDirection(): S[abs(S) < abs(t(S))] <- 0
        S = A0 .* ~(abs(A0) < abs(A0.'));
    else
        S = A0 .* (abs(A0) > abs(A0.'));
    end
    A0 = (1 - lambdav) * A0 + lambdav * S;
end
A0 = A0 - diag(diag(A0));
% Row k now lists the targets of gene k
A0 = A0.';

nOutgoing = nnz(A0(idx, :) ~= 0);
if reference
    if nOutgoing == 0
        warning('ten:i_knk:noOutgoingEdges', ...
            ['%s has no outgoing edges in the WT network; the knockout does ' ...
             'not change the network and the differential regulation results ' ...
             'reflect only numerical noise.'], genelist(idx));
    end
elseif nOutgoing < 10
    warning('KO gene (%s) has no link or too few links with other genes.', ...
        genelist(idx));
    T = table();
    return;
end

A1 = A0;
A1(idx, :) = 0;

if reference
    [aln0, aln1] = ten.i_ma(A0, A1, 2, "unnormalized", false);
    T = ten.i_dr(aln0, aln1, genelist, dosort, true);
    return
end

a = A0(idx, :);
b = A1(idx, :);
yn = abs(a-b) <= eps(max(abs(a), abs(b)));
yn = all(yn); % scalar logical output

if yn
    drdist = nan(length(genelist), 1);
    T = table(drdist);
else
    % Normalized Laplacian (the ten.i_ma default); see its header for the
    % Knk benchmark and when "unnormalized" is the better choice.
    [aln0, aln1] = ten.i_ma(A0, A1, 2);
    T = ten.i_dr(aln0, aln1, genelist, dosort);
end
end
