function q = e_fdr(p)
%E_FDR  Benjamini-Hochberg adjusted p-values. NaN-safe, no toolbox needed.
%
%   q = PKG.E_FDR(p) returns the Benjamini-Hochberg adjusted p-value for
%   every entry of P, preserving its shape. Entries that are not finite --
%   a gene that could not be tested at all -- come back as NaN and are
%   excluded from the family size, because a hypothesis that was never
%   tested is not one of the hypotheses being corrected over.
%
%   WHY THIS EXISTS. Five functions in this toolbox carried the same block:
%
%       if exist('mafdr.m', 'file')
%           p_val_adj = mafdr(p_val, 'BHFDR', true);
%       else
%           [~, ~, ~, p_val_adj] = pkg.e_fdr_bh(p_val);
%       end
%
%   The two branches are the same algorithm -- on p-values with no NaN they
%   agree to 2e-16 -- but they part company on NaN, which is exactly what a
%   rank-sum test returns for a gene with no counts in either group, and
%   what a per-cell-type differential expression run produces in bulk.
%   MAFDR drops those from the family; PKG.E_FDR_BH keeps them. So the same
%   GUI command gave a different answer depending on whether the
%   Bioinformatics Toolbox happened to be installed: on 20000 genes with
%   8000 undetected, the correction ran over 20000 hypotheses instead of
%   12000 and inflated every adjusted p-value by a factor of 1.67. The
%   no-toolbox path was the conservative one, so genes that should have
%   been called were not.
%
%   PKG.E_FDR_BH also returns values above 1 in that case -- 1796 of them
%   in the run above, up to 1.67. Its own help says adjusted p-values can
%   exceed 1, which is true of the raw step-up expression but is not a
%   probability, and every caller here compares the result against 0.05.
%
%   This function is the single behaviour: it reproduces MAFDR(...,
%   'BHFDR', true) exactly, on NaN-free and NaN-carrying input alike, and
%   needs no toolbox to do it.
%
%   PKG.E_FDR_BH is left in place but nothing in the toolbox calls it any
%   more. Every one of its eight call sites took only its fourth output,
%   the adjusted p-value; not one used the rejection vector H or the
%   critical p-value that are its other reason to exist. Two of those sites
%   carried comments warning that it "returns h FIRST and adjusted p
%   FOURTH", and a pipeline note spells out that taking output 1 as a
%   q-value "silently inverts the test and retains almost everything --
%   which looks like a working filter, because a number did come out".
%   This function has one output, so that trap cannot be sprung.
%
%   INPUT:
%     p - p-values, any shape. Non-finite entries are treated as untested.
%
%   OUTPUT:
%     q - adjusted p-values in [0, 1], NaN where P was, same shape as P.
%
%   See also PKG.E_FDR_BH, MAFDR, SC_DEG, SC_DPG, SC_PICKMARKERS.

arguments
    p {mustBeNumeric}
end

q = NaN(size(p));
tested = isfinite(p);
numTested = nnz(tested);
if numTested == 0
    return;
end

pv = double(p(tested));
[sorted, order] = sort(pv(:));

% Step-up: q_(i) = min over j >= i of (m/j)*p_(j), then capped at 1. The
% reverse cumulative minimum is what enforces monotonicity, without which
% a large p-value early in the sequence could be reported as more
% significant than a smaller one after it.
adjusted = min(1, cummin(sorted.*(numTested./(1:numTested).'), "reverse"));

restored = zeros(numTested, 1);
restored(order) = adjusted;
q(tested) = restored;

end
