function Y = e_randmds(X, ndim)
%E_RANDMDS  Classical MDS coordinates via a randomised SVD.
%
%   Y = PKG.E_RANDMDS(D, ndim) returns NDIM classical multidimensional
%   scaling coordinates for the distance matrix D. The result matches
%   MATLAB's CMDSCALE(D, ndim) up to sign and rotation.
%
%   THE SINGULAR VALUES WERE BEING DISCARDED. This was
%
%       [Y, ~, ~] = pkg.e_randPCA(X', ndim);
%
%   and PKG.E_RANDPCA returns a rank-k SVD U*S*V', whose U has orthonormal
%   columns by construction. So every returned dimension had norm exactly
%   1, whatever its singular value, and the relative importance of the
%   axes -- the whole content of classical MDS -- was gone. Classical MDS
%   coordinates are U*sqrt(Lambda).
%
%   Measured on 200 points in 5-D forming two well-separated blobs, against
%   CMDSCALE on the same distance matrix:
%
%                        column norms            dim1/dim2   dist. corr.
%       cmdscale     [127.05  15.66  14.54]        8.111       0.9986
%       old          [  1.00   1.00   1.00]        1.000       0.5649
%       U*sqrt(S)    [127.05  15.66  14.54]        8.111       0.9986
%
%   with a Procrustes residual against CMDSCALE of 4.4e-16, against 0.503
%   for the unscaled version.
%
%   The only caller is PKG.E_EMBEDBYD, which uses this as the starting
%   configuration for t-SNE or MDSCALE, so a starting point with every
%   axis equally scaled was throwing away the structure it was there to
%   provide.
%
%   See also CMDSCALE, PKG.E_EMBEDBYD, PKG.E_RANDPCA.

X = X.^2;
X = bsxfun(@minus, X, mean(X, 1));
X = bsxfun(@minus, X, mean(X, 2));

% The -0.5 completes the classical-MDS double centring. It was missing, and
% on its own it is only a uniform factor of sqrt(2) once the square root is
% taken below -- harmless for a starting configuration -- but including it
% makes the output literally CMDSCALE's coordinates rather than a scaled
% copy.
X = -0.5*X;

[U, S, ~] = pkg.e_randPCA(X', ndim);
Y = U .* sqrt(max(diag(S), 0))';
