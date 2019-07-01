function [score]=fun_cmp_cluster(true_labels,cluster_labels,varargin)

% https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bby076/5077112

p = inputParser;
defaultType = 'nmi';
validTypes = {'nmi','ari','fm','jaccard'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'true_labels',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,true_labels,varargin{:})

addRequired(p,'cluster_labels',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,cluster_labels,varargin{:})


pw1=fileparts(which(mfilename));
switch p.Results.type
    case 'nmi'        
        pth=fullfile(pw1,'thirdparty/SIMLR');
        addpath(pth);
        pth=fullfile(pw1,'thirdparty/SIMLR/src');
        addpath(pth);
        score = Cal_NMI(true_labels, cluster_labels);
        fprintf('The NMI value is %f\n', NMI_i);
    case 'ari'
        % ----
    case 'fm'
        % ----
    case 'jaccard'
        % ----

end
