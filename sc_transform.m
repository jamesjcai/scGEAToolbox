function [X]=sc_transform(X,varargin)

p = inputParser;
defaultType = 'pearsonresiduals';
validTypes = {'pearsonresiduals','csndm','bigscale','sct'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,varargin{:})

   
switch p.Results.type
    case 'pearsonresiduals'
        % analytic Pearson residuals
        % https://doi.org/10.1101/2020.12.01.405886
        % https://gist.github.com/hypercompetent/51a3c428745e1c06d826d76c3671797c
        
        u=(sum(X,2)*sum(X,1))./sum(X(:));
        s=sqrt(u+(u.^2)./100);
        X=(X-u)./s;
        n=size(X,2);
        % clip to sqrt(n);
        sn=sqrt(n);
        X(X>sn)=sn;
        X(X<-sn)=-sn;        
    case 'glmpca'
        
    case 'normalization_sqrt'
        % https://twitter.com/hippopedoid/status/1337028817219620864?s=20        
    case 'Hafemeister & Satija 2019'
        % https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
    case 'csndm'
        [X]=run.CSN_transform(X);
%     case 'bigscale'
%         pw1=fileparts(which(mfilename));
%         pth=fullfile(pw1,'thirdparty/bigSCale');
%         addpath(pth);
%         % model=1. Log(x), then each row (gene) normalized between [-5:5]
%         [X]=transform_bigscale(X);
    case 'sct'
        % sc_sct
end
