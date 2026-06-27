function [Y] = e_tsnebyd(D, ndim, InitialY)
if nargin < 2, ndim = 2; end
if nargin < 3, InitialY = []; end

% tsne has no 'precomputed' Distance option, so feed it row indices and use
% a closure over D to look up the precomputed pairwise distance. ZI is a
% 1x1 scalar (the index i), ZJ is m x 1 (indices [j1; ...; jm]); D(ZJ, ZI)
% returns an m x 1 column, which is what tsne expects.
N = size(D, 1);
X = (1:N)';
distFun = @(ZI, ZJ) D(ZJ, ZI);

if ~isempty(InitialY)
    Y = tsne(X, 'Distance', distFun, ...
        'NumDimensions', ndim, 'InitialY', InitialY);
else
    Y = tsne(X, 'Distance', distFun, ...
        'NumDimensions', ndim);
end
end
