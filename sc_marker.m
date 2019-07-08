function [idx]=sc_marker(X,genelist,c,varargin)

p = inputParser;
defaultType = 'soptsc';
validTypes = {'soptsc'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addRequired(p,'genelist',@isstring);
addRequired(p,'c',@isnumeric);
addOptional(p,'type',defaultType,checkType)
addOptional(p,'numofmarkers',10,@(x) (x > 0) && isinteger(x))
addOptional(p,'plotit',false,@islogical)
parse(p,X,genelist,c,varargin{:})

switch p.Results.type
    case 'soptsc'
        idx=GC_htmp_DE(X,genelist,c,...
            p.Results.numofmarkers,p.Results.plotit);
end
