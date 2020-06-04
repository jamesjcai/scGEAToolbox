function [A]=sc_grnetwork(X,varargin)

p = inputParser;
defaultType = 'pcnetpar';
validTypes = {'pcnetpar','pcnet','genie3'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,varargin{:})

   
switch p.Results.type
    case 'pcnetpar'
        [A]=sc_pcnetpar(X);
    case 'pcnet'
        [A]=sc_pcnet(X);
    case 'genie3'
        [A]=run_genie3(X,[],true);
end
