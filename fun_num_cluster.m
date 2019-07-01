function [nC]=fun_num_cluster(X,varargin)

p = inputParser;
defaultType = 'simlr';
validTypes = {'simlr','soptsc'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,varargin{:})

pw1=fileparts(which(mfilename));
switch p.Results.type
    case 'simlr'        
        pth=fullfile(pw1,'thirdparty/SIMLR');
        addpath(pth);
        pth=fullfile(pw1,'thirdparty/SIMLR/src');
        addpath(pth);
        [K1, K2] = Estimate_Number_of_Clusters_SIMLR(X',2:10);
        [~,i]=min(K2);
        nC=i+1;
    case 'soptsc'        
        pth=fullfile(pw1,'thirdparty/SoptSC');
        addpath(pth);
        pth=fullfile(pw1,'thirdparty/SoptSC/NNDSVD');
        addpath(pth);
        pth=fullfile(pw1,'thirdparty/SoptSC/symnmf2');
        addpath(pth);
        
        
        
        realdata = X;
        realdata = realdata-min(realdata(:));
        realdata = realdata./max(realdata(:));

        [~,n] = size(realdata);
        for i = 1:n
            realdata(:,i) = realdata(:,i)/norm(realdata(:,i));
        end
        lambda = 0.5;
        W = SimilarityM(realdata,lambda,X);
        WB = W;
        n = size(W,1);
        D = diag(WB*ones(n,1));
        Prw = eye(size(W)) - D^(-1/2)*WB*D^(-1/2);
        if n>=1000
            No_eigs = 100;
            all_eigs = real(eigs(Prw,No_eigs,'sm'));
        else
            all_eigs = real(eig(Prw));
        end

        ZZ = sort(abs(real(all_eigs)));        
        No_cluster1 = length(find(ZZ<=0.01));

        % Determinning the number of clusters
        eigenvalues = [];
        if isempty(nC)
            [eigenvalues,No_cluster] = Num_cluster(W,No_cluster1);
            nC = No_cluster;
        end
        
end
