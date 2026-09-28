function [v, mu] = e_rowvar(X, nanflag)
%E_ROWVAR Variance and mean of each row, fast on sparse matrices.
%   [V, MU] = pkg.e_rowvar(X) returns var(X, 0, 2) and mean(X, 2) as full
%   column vectors. [V, MU] = pkg.e_rowvar(X, "omitnan") ignores NaNs, as
%   var(X, 0, 2, "omitnan") does. Take SQRT(V) for the standard deviation.
%
%   For a dense X this is VAR and MEAN themselves. For a sparse X, VAR
%   works through every zero; this takes the row sums and sums of squares
%   with sparse arithmetic, which touches only the stored entries. At
%   20000 genes x 30000 cells std(X, 0, 2) took 7.0 s.
%
%   Sums of squares alone (E[x^2] - E[x]^2) lose the variance to
%   cancellation when a row's mean is large next to its spread, and can
%   leave a small positive number on a constant row where VAR returns
%   exactly 0 -- which a caller taking log(cv) would carry through. So any
%   row where the squared deviations come to under 0.1% of the sum of
%   squares, and any row holding a NaN, is recomputed with VAR's own
%   two-pass formula (the mean first, then the squared deviations from it)
%   over its stored entries: the zeros all deviate by the same -MU, so they
%   add (number of zeros) * MU^2. The result agrees with VAR to about 1e-13
%   relative.
%
%   See also var, std, mean, sc_analyticfit, sc_hvg.

arguments
    X {mustBeNumeric}
    nanflag (1,1) string {mustBeMember(nanflag, ["includenan", "omitnan"])} = "includenan"
end

if ~issparse(X)
    v = var(X, 0, 2, nanflag);
    mu = mean(X, 2, nanflag);
    return;
end

N = size(X, 2);
omit = nanflag == "omitnan";

% Fast path: sums and sums of squares by sparse arithmetic.
sx = full(sum(X, 2));
sq = full(sum(X.^2, 2));
mu = sx/N;
ss = sq - sx.*mu;
v = max(ss, 0)/max(N - 1, 1);

% Rows the fast path may get wrong: cancellation, or a NaN in the row.
redo = ~(ss > 1e-3*sq) & sq > 0;
redo = redo | isnan(sx) | isnan(sq);
if any(redo)
    rows = find(redo);
    [vr, mur] = i_twopass(X(rows, :), omit);
    v(rows) = vr;
    mu(rows) = mur;
end
end


function [v, mu] = i_twopass(X, omit)
% VAR's two-pass formula over the stored entries of sparse X.
[G, N] = size(X);

% Pass 1: per-row sums, nonzero counts and NaN counts.
s = zeros(G, 1);
nz = zeros(G, 1);
nn = zeros(G, 1);
blocks = i_blocks(X);
for b = 1:size(blocks, 1)
    [i, ~, x] = find(X(:, blocks(b, 1):blocks(b, 2)));
    i = i(:);   % FIND gives rows for a one-row X
    x = double(x(:));
    bad = isnan(x);
    nz = nz + accumarray(i, 1, [G, 1]);
    nn = nn + accumarray(i, bad, [G, 1]);
    if omit
        s = s + accumarray(i(~bad), x(~bad), [G, 1]);
    else
        s = s + accumarray(i, x, [G, 1]);
    end
end

n = N*ones(G, 1);
if omit
    n = n - nn;                         % NaNs are stored entries, not zeros
end
mu = s./n;
n0 = N - nz;                            % the zeros

% Pass 2: squared deviations from the mean, nonzeros plus the zero block.
ss = n0.*mu.^2;
for b = 1:size(blocks, 1)
    [i, ~, x] = find(X(:, blocks(b, 1):blocks(b, 2)));
    i = i(:);   % FIND gives rows for a one-row X
    x = double(x(:));
    if omit
        ok = ~isnan(x);
        i = i(ok);
        x = x(ok);
    end
    ss = ss + accumarray(i, (x - mu(i)).^2, [G, 1]);
end

% VAR divides by n - 1, and returns 0 for a single observation.
v = ss./max(n - 1, 1);
v(n == 0) = NaN;
if ~omit
    v(nn > 0) = NaN;
    mu(nn > 0) = NaN;
end
end


function blocks = i_blocks(X)
% Column ranges holding about 1e7 stored entries each, to bound the memory
% of FIND's three output vectors.
N = size(X, 2);
perCol = max(1, nnz(X)/max(N, 1));
step = max(1, floor(1e7/perCol));
starts = (1:step:N).';
blocks = [starts, min(starts + step - 1, N)];
if isempty(blocks)
    blocks = zeros(0, 2);
end
end
