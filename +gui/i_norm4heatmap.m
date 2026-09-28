function [Y] = i_norm4heatmap(Y, dim, methodid, doclip)

% Y - gene-by-cell matrix
% dim = 2 - by row
if nargin < 4, doclip = true; end
if nargin < 3 || isempty(methodid), methodid = 1; end
if nargin < 2 || isempty(dim), dim = 2; end

% One method at a time. A multi-select dialog upstream used to hand this a
% two-element METHODID, and all the SWITCH below could say about it was
% "SWITCH expression must be a scalar or a character vector", which names
% neither the dialog nor the argument.
if ~isscalar(methodid)
    error('gui:i_norm4heatmap:notScalar', ...
        ['METHODID must be a single method, not %d of them. ', ...
        'Select one entry in the normalization method dialog.'], ...
        numel(methodid));
end


switch methodid
    case 1
        Y = normalize(Y, dim, 'zscore');
    case 2
        Y = normalize(Y, dim, 'zscore', 'robust');
    case 3
        Y = normalize(Y, dim, 'norm', 2);
    case 4
        Y = normalize(Y, dim, 'norm', Inf);
    case 5
        Y = normalize(Y, dim, 'scale', 'std');
    case 6
        Y = normalize(Y, dim, 'scale', 'mad');
    case 7
        Y = normalize(Y, dim, 'scale', 'first');
    case 8
        Y = normalize(Y, dim, 'scale', 'iqr');
    case 9
        Y = normalize(Y, dim, 'range', [0 1]);
    case 10
        Y = normalize(Y, dim, 'range', [1 10]);
    case 11
        Y = normalize(Y, dim, 'center', 'mean');
    case 12
        Y = normalize(Y, dim, 'center', 'median');
    case 13
        Y = normalize(Y, dim, 'center', 'scale');
    case 14
        Y = normalize(Y, dim, 'medianiqr');
    otherwise
        error('gui:i_norm4heatmap:badMethod', ...
            ['Unknown normalization method %s. ', ...
            'Expected an integer from 1 to 14.'], mat2str(methodid));
end


if doclip
    q = quantile(Y(:), [0.05, 0.95]);
    Y(Y < q(1)) = q(1);
    Y(Y > q(2)) = q(2);
end

end
