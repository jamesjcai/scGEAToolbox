function [X]=sc_impute(X,varargin)
% Imputation 
% 
% See alos: SC_TRANSFORM

p = inputParser;
defaultType = 'MAGIC';
validTypes = {'MAGIC','McImpute','SAVER'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,varargin{:})

   
switch upper(p.Results.type)
    case 'MAGIC'
        [X]=run.MAGIC(X,true);
    case 'MCIMPUTE'
        [X]=run.McImpute(X,true);
    case 'SAVER'
        [X]=run.SAVER(X);
end
end