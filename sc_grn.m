function [A]=sc_grn(X,varargin)
% construct single-cell gene regulatory network (scGRN)

p = inputParser;
defaultType = 'pcnetpar';
validTypes = {'pcnet','pcnetpar','pcnetdenoised','genie3'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,varargin{:})

   
switch p.Results.type
    case 'pcnet'
        [A]=sc_pcnet(X);
    case 'pcnetpar'
        [A]=sc_pcnetpar(X);
    case 'pcnetdenoised'
        [A]=sc_pcnetdenoised(X);
    case 'genie3'
        [A]=run.GENIE3(X,[],true);
end
end