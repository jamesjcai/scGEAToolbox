function [X]=sc_norm(X,varargin)

p = inputParser;
defaultType = 'libsize';
validTypes = {'libsize','deseq'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,varargin{:})

if (max(sum(X))-min(sum(X)))<0.1
    warning('Input X may have already been normalized.');
end
    
switch p.Results.type
    case 'libsize'
        [X]=norm_libsize(X);
    case 'deseq'
        [X]=norm_deseq(X);
end
