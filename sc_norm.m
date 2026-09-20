function [X] = sc_norm(X, varargin)

p = inputParser;
defaultType = 'libsize';
validTypes = {'libsize', 'deseq', 'shiftedclr'};
checkType = @(x) any(validatestring(x, validTypes));

addRequired(p, 'X', @isnumeric);
addOptional(p, 'type', defaultType, checkType);
parse(p, X, varargin{:});

% Near-identical library sizes across cells suggest X was normalized
% already. Guard the degenerate cases that also give a zero spread but say
% nothing about normalization: a single cell, and an all-zero submatrix
% (e.g. a gene subset with no counts in the selected cells).
sx = sum(X, 1, 'omitnan');
if numel(sx) > 1 && max(sx) > 0 && (max(sx)-min(sx))/max(sx) < 1e-6
    warning('sc_norm:AlreadyNormalized', ...
        'Input X may have already been normalized.');
end

switch p.Results.type
    case 'libsize'
        [X] = pkg.norm_libsize(X);
    case 'deseq'
        [X] = pkg.norm_deseq(X);
    case 'shiftedclr'
        [X] = pkg.norm_shiftedclr_scclr(X);
    otherwise
        error('sc_norm:InvalidType', 'Unknown normalization type: %s', p.Results.type);
end
end
