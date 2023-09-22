function [Y] = e_tsnebyd(D, ndim, InitialY)
if nargin < 2, ndim = 2; end
if nargin < 3, InitialY = []; end
if ~isempty(InitialY)
    Y = tsne(D, 'Distance', @ExampleDistFunc, ...
        'NumDimensions', ndim, 'InitialY', InitialY);
else
    Y = tsne(D, 'Distance', @ExampleDistFunc, ...
        'NumDimensions', ndim);
end
end

function D2 = ExampleDistFunc(Z1, ZJ)
n = min([50, length(Z1)]);
z1 = Z1(1:n);
zj = ZJ(:, 1:n);
[~, index] = ismember(z1, zj, 'rows');
D2 = ZJ(:, index);
end
