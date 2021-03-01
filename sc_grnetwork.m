function [A]=sc_grnetwork(X,varargin)

p = inputParser;
defaultType = 'pcnetpar';
validTypes = {'pcnet','pcnetpar','genie3'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,varargin{:})

   
switch p.Results.type
    case 'pcnet'
        [A]=sc_pcnetpar(X);
    case 'pcnetpar'
        [A]=sc_pcnet(X);
    case 'genie3'
        [A]=run.GENIE3(X,[],true);
end
