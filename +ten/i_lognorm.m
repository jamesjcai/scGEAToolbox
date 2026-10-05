function X = i_lognorm(X)
%I_LOGNORM  Median library-size normalisation followed by log1p.
%   X = TEN.I_LOGNORM(X) scales each column (cell) of the genes-by-cells
%   count matrix X to the median library size and applies log1p. Empty cells
%   are left at zero rather than divided by zero.
%
%   The one normalisation for the scTenifoldXct code: the native cores and
%   xctmain functions apply it before building their networks, TEN.I_XCTW12
%   before computing its metrics, and the +run/py_scTenifoldXct wrappers
%   before building the networks they hand to Python.
%
%   It differs from log1p(sc_norm(X)). sc_norm scales every cell to 1e4, so
%   on shallow data (median library size far below 1e4) the zero/non-zero
%   gap dominates after log1p and per-cell depth comes back as
%   co-expression; scaling to the median keeps each cell near its own depth.
%   sc_norm also turns empty cells into NaN.
%
% see also: TEN.I_PCNET, TEN.I_XCTCORE, TEN.I_XCTW12

X = double(X);
colSums = sum(X, 1);
colSums(colSums == 0) = 1;    % avoid /0 for empty cells
X = X./colSums.*median(colSums);
X = log1p(X);
end
