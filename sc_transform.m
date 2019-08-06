function [X]=sc_transform(X,varargin)

p = inputParser;
defaultType = 'csndm';
validTypes = {'csndm','bigscale'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,varargin{:})

   
switch p.Results.type
    case 'csndm'
        [X]=run_csndm(X);
    case 'bigscale'
        pw1=fileparts(which(mfilename));
        pth=fullfile(pw1,'thirdparty/bigSCale');
        addpath(pth);
        % model=1. Log(x), then each row (gene) normalized between [-5:5]
        [X]=transform_bigscale(X);
end
