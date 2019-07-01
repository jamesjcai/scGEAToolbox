function [C]=sc_cluster(X,varargin)

p = inputParser;
defaultType = 'simlr';
validTypes = {'simlr','simlr_pearson','soptsc','alona'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,varargin{:})

switch p.Results.type
    case 'simlr'
        [C]=run_simlr(X);
    case 'simlr_pearson'
        [C]=run_simlr(X);
    case 'soptsc'
        [C]=run_soptsc(X);
    case 'alona'
        [C]=sc_alona(X);
end
