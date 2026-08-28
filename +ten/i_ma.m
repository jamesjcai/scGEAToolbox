function [aln0, aln1] = i_ma(A0, A1, ndim)
% MA - manifold alignment

if nargin < 3, ndim = 30; end
mu = 0.9;

W1 = A0 + 1;
W2 = A1 + 1;

W12 = eye(size(W1, 2), size(W2, 2));
mu = mu * (sum(W1(:)) + sum(W2(:))) / (2 * sum(W12(:)));
W = [W1, mu * W12; mu * W12', W2];

D = sum(abs(W));
L = diag(D) - W;
L = 0.5 * (L + L.');    % symmetric by construction; make it exactly so

V = i_lomodes(L, ndim);

p1 = size(W1, 1);
aln0 = V(1:p1, :);
aln1 = V(p1+1:end, :);
end


function V = i_lomodes(L, ndim)
% The ndim eigenvectors of smallest eigenvalue, skipping the null space.
%
% EIGS IS ITERATIVE AND DOES NOT ALWAYS CONVERGE ON THIS LAPLACIAN. When it
% does not it warns, returns NaN for the eigenvalues that failed, and the
% columns left after the null-space filter can fall below ndim - which used
% to surface as "Index in position 2 exceeds array bounds" from the line
% below, several hours into a scTenifoldNet run. Two things also make a
% partial result unusable rather than merely short: the converged subset is
% not guaranteed to be the smallest ones, and a disconnected aligned graph
% puts more than ndim eigenvalues at zero, so ndim*2 is not always enough
% cushion.
%
% So: try eigs, and on anything less than full convergence redo it with a
% dense symmetric eig, which is exact and affordable because L is dense
% already. That costs seconds at the usual sizes and grows as n^3, so it is
% a fallback rather than the default.
[V, d] = i_eigs(L, ndim);
if isempty(V) || any(isnan(d)) || sum(d >= 1e-8) < ndim
    if isempty(V)
        why = 'eigs failed';
    elseif any(isnan(d))
        why = sprintf('%d of %d eigenvalues did not converge', ...
            sum(isnan(d)), numel(d));
    else
        why = sprintf('only %d eigenvalue(s) above the null-space cutoff', ...
            sum(d >= 1e-8));
    end
    warning('ten:i_ma:denseFallback', ...
        '%s; recomputing the %d smallest modes with a dense eig (n = %d).', ...
        why, ndim, size(L, 1));
    [V, d] = eig(full(L), 'vector');
    [d, ind] = sort(d);
    V = V(:, ind);
end

V = V(:, d >= 1e-8);
assert(size(V, 2) >= ndim, 'ten:i_ma:tooFewModes', ...
    ['only %d eigenvalue(s) above 1e-8 but %d modes needed - the aligned ' ...
     'graph is more disconnected than the alignment assumes.'], ...
    size(V, 2), ndim);
V = V(:, 1:ndim);
end


function [V, d] = i_eigs(L, ndim)
% ndim*2 leaves room for the null space. Returns empty V if eigs errors
% outright, so the caller can fall back on that too.
V = []; d = [];
try
    [V, D] = eigs(L, ndim*2, 'smallestreal');
catch
    return
end
d = diag(D);
[d, ind] = sort(d);
V = V(:, ind);
end
