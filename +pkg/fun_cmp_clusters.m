function [score]=fun_cmp_clusters(true_labels,cluster_labels,varargin)
%Compare two clusters
% https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bby076/5077112

p = inputParser;
defaultType = 'nmi';
validTypes = {'nmi','ari','fm','jaccard'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'true_labels',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,true_labels,varargin{:})


pw1=fileparts(mfilename('fullpath'));
switch p.Results.type
    case 'nmi'        
        pth=fullfile(pw1,'thirdparty','SIMLR');
        if ~(ismcc || isdeployed), addpath(pth); end
        pth=fullfile(pw1,'thirdparty','SIMLR','src');
        if ~(ismcc || isdeployed), addpath(pth); end
        score = Cal_NMI(true_labels, cluster_labels);
        fprintf('The NMI value is %f\n', score);
    case 'ari'
        % ----
        pth=fullfile(pw1,'+run','thirdparty','SinNLRR');
        if ~(ismcc || isdeployed), addpath(pth); end
        [AR,~,~,~]=Cal_ARI(true_labels, cluster_labels);
        score=AR;
        fprintf('The ARI value is %f\n', score);
        
    case 'fm'
        % ----
    case 'jaccard'
        % ----

end
end