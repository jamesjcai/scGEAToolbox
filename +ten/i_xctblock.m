function [W, info] = i_xctblock(W11, W22, W12, mu, offset)
%I_XCTBLOCK  Assemble the two-block alignment graph, with the GRN offset.
%   W = TEN.I_XCTBLOCK(W11, W22, W12, MU, OFFSET) builds
%
%       W = [ W11+OFFSET,        s*W12 ]
%           [ s*W12',       W22+OFFSET ],   s = MU*(sum|W11'|+sum|W22'|)
%                                               /(2*sum|W12|)
%
%   where W11' and W22' are the offset blocks, so the scale factor is computed
%   from the offset networks exactly as core.py:432 does.
%
%   THE OFFSET IS ADDED TO EVERY ELEMENT, not to the diagonal. This is the
%   single most-confused line in this codebase's history: the reference and
%   ten.i_ma both write A+1, ten.xct.xctmain rendered it as W11+speye(ng), and
%   ten.sctenifoldxct dropped it while keeping a comment claiming it was there.
%   A diagonal shift is a no-op for the Laplacian - adding c to W_ii adds c to
%   both d_i and W_ii, and L_ii = d_i - W_ii is unchanged - so the two are not
%   interchangeable. Keeping the assembly in one function is what stops that
%   divergence recurring.
%
%   WHAT THE OFFSET DOES. It makes every block dense and non-negative, which
%   (a) connects the graph, collapsing the component count to 1, (b) makes the
%   Laplacian positive semidefinite, and (c) makes the sum(W) and sum(abs(W))
%   degree conventions identical, so MATLAB and the reference agree. Measured
%   on the standard fixture: 168 components and 43 negative eigenvalues without
%   it, 1 component and none with it.
%
%   It is NOT a uniform spectral shift. Adding c to every element of the whole
%   matrix would give L' = L + c*(n*I - J), which preserves eigenvectors and
%   would merely relocate the degenerate directions rather than remove them.
%   That is not what happens here, because the offset applies only to the two
%   diagonal blocks and because the scale factor s is recomputed from the
%   offset networks, which changes the cross-block coupling by orders of
%   magnitude. Measured: the offset spectrum's bottom non-trivial modes retain
%   0.2% of their energy in the un-offset null space.
%
%   INPUTS:
%     W11, W22 - NG-by-NG symmetric within-type GRN blocks
%     W12      - NG-by-NG cross-type correspondence block
%     MU       - cross-type coupling weight
%     OFFSET   - constant added to every element of W11 and W22; 0 disables
%
%   OUTPUTS:
%     W    - 2*NG-by-2*NG block matrix, dense when OFFSET is nonzero
%     info - struct with .offset, .mu_scale, .dense
%
%   MEMORY. A nonzero offset makes W dense: 2*NG-by-2*NG doubles, which is
%   ~110 MB at NG=1856 and ~3.2 GB at NG=10000. The reference has the same
%   ceiling and relies on the caller subsetting genes; so does this.
%
% see also: TEN.I_XCTCORE, TEN.I_XCTW12, TEN.XCT.XCTMAIN

arguments
    W11 {mustBeNumeric}
    W22 {mustBeNumeric}
    W12 {mustBeNumeric}
    mu (1, 1) double
    offset (1, 1) double = 1
end

ng = size(W11, 1);

if offset ~= 0
    n = 2*ng;
    if n > 30000
        warning("TEN:I_XCTBLOCK:LargeDense", ...
            "A nonzero GRN offset makes the %d-by-%d alignment graph dense " + ...
            "(~%.1f GB). Subset genes first, or pass offset=0.", ...
            n, n, (n^2*8)/2^30);
    end
    W11 = full(double(W11)) + offset;
    W22 = full(double(W22)) + offset;
    W12 = full(double(W12));
    isDense = true;
else
    W11 = sparse(W11);
    W22 = sparse(W22);
    isDense = false;
end

% Scale factor from the offset networks, as core.py:432. With a nonzero offset
% every block is non-negative, so these abs sums equal the reference's plain
% sums; with offset 0 the abs form preserves the previous MATLAB behaviour on a
% signed GRN.
w12_sum = full(sum(abs(W12(:))));
if w12_sum == 0
    mu_scale = 0;
else
    mu_scale = mu*(full(sum(abs(W11(:)))) + full(sum(abs(W22(:)))))/(2*w12_sum);
end

W = [W11,               mu_scale.*W12; ...
     mu_scale.*W12',    W22];

info = struct("offset", offset, "mu_scale", mu_scale, "dense", isDense);

end
