function [X]=sc_transform(X,varargin)

% https://www.biorxiv.org/content/10.1101/2021.06.24.449781v1.full
% acosh transformation based on the delta method
% shifted logarithm (log(x + c)) with a pseudo-count c, so that it approximates the acosh transformation
% randomized quantile and Pearson residuals

p = inputParser;
defaultType = 'PearsonResiduals';
validTypes = {'PearsonResiduals','kNNSmoothing','FreemanTukey','csndm','sct'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,varargin{:})

   
switch lower(p.Results.type)
    case 'pearsonresiduals'
        % analytic Pearson residuals
        % https://doi.org/10.1101/2020.12.01.405886
        % https://gist.github.com/hypercompetent/51a3c428745e1c06d826d76c3671797c
        
        u=(sum(X,2)*sum(X,1))./sum(X(:));
        s=sqrt(u+(u.^2)./100);
        X=(X-u)./s;
        X(isnan(X))=0;
        n=size(X,2);
        % clip to sqrt(n);
        sn=sqrt(n);
        X(X>sn)=sn;
        X(X<-sn)=-sn;       
    case 'knnsmoothing'
        % K-nearest neighbor smoothing for high-throughput single-cell RNA-Seq data
        % https://doi.org/10.1101/217737
        % 
        X=knn_smooth(X,5,10);
        
    case 'glmpca'
        
    case 'normalization_sqrt'
        % https://twitter.com/hippopedoid/status/1337028817219620864?s=20        
    case 'Hafemeister & Satija 2019'
        % https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
    case 'csndm'
        [X]=run.csndm_trans(X);
%     case 'bigscale'
%         pth=fullfile(pw1,'thirdparty/bigSCale');
%         addpath(pth);
%         % model=1. Log(x), then each row (gene) normalized between [-5:5]
%         [X]=transform_bigscale(X);
    case 'sct'
        % sc_sct
    case 'freemantukey'
        % https://github.com/flo-compbio/monet/blob/master/monet/util/expression.py
        % Applies the Freeman-Tukey transformation to stabilize variance."
        % https://www.biorxiv.org/content/10.1101/2020.06.08.140673v2.full
        % https://www.nature.com/articles/nmeth.2930
        X=sc_norm(X,'type','deseq');
        X=sqrt(X)+sqrt(X+1);
end
end