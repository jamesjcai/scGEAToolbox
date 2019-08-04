function [X]=sc_transform(X,varargin)

p = inputParser;
defaultType = 'csndm';
validTypes = {'csndm'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,varargin{:})

   
switch p.Results.type
    case 'csndm'
        [X]=run_csndm(X);
end
