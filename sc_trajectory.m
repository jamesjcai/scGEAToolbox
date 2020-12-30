function [t]=sc_trajectory(X,varargin)

p = inputParser;
defaultType = 'splinefit';
validTypes = {'splinefit','tscan'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
addOptional(p,'plotit',false,@islogical);
parse(p,X,varargin{:})
plotit=p.Results.plotit;

switch p.Results.type
    case 'splinefit'
        s=run.PHATE(X,3,plotit,false);
        [t]=i_pseudotime_by_splinefit(s,1,plotit);
    case 'tscan'
        t=sc_tscan(X,'plotit',true);        
end
if size(t,2)~=1, t=t'; end
end

