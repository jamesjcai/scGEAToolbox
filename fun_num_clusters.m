function [optimk]=fun_num_clusters(X,varargin)
%Estimate number of clusters
p = inputParser;
defaultType = 'simlr';
validTypes = {'simlr','soptsc','sc3'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,varargin{:})


pw1=fileparts(mfilename('fullpath'));
switch p.Results.type
    case 'simlr'        
        pth=fullfile(pw1,'+run','thirdparty','SIMLR');
        addpath(pth);
        pth=fullfile(pw1,'+run','thirdparty','SIMLR','src');
        addpath(pth);
        [K1, K2] = Estimate_Number_of_Clusters_SIMLR(X',2:10);
        [~,i]=min(K2);
        optimk=i+1;
    case 'soptsc'        
        pth=fullfile(pw1,'+run','thirdparty','SoptSC');
        addpath(pth);
        pth=fullfile(pw1,'+run','thirdparty','SoptSC','NNDSVD');
        addpath(pth);
        pth=fullfile(pw1,'+run','thirdparty','SoptSC','symnmf2');
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
        %if isempty(optimk)
            [eigenvalues,No_cluster] = Num_cluster(W,No_cluster1);
            optimk = No_cluster;
        %end
    case 'sc3'
        %% estimate k
        % X=log2(X+1);
        Dis=squareform(pdist(X'));
        A=exp(-Dis./max(Dis(:)));   % adjacency matrix
        xD=diag(sum(A).^-0.5);  % D=diag(sum(A)); % d(i) the degree of node i
        xA=xD*A*xD;             % normalized adjacenty matrix
        L=eye(size(A,1))-xA;    % also L=xD*(D-A)*xD 

        % see https://people.orie.cornell.edu/dpw/orie6334/lecture7.pdf
        % see https://en.wikipedia.org/wiki/Laplacian_matrix#Symmetric_normalized_Laplacian_2

        [V,D]=eig(L);
        [~,ind] = sort(diag(D));
        Ds = D(ind,ind);
        Vs = V(:,ind);

        clust=zeros(size(Vs,1),6);
        for i=1:6
            clust(:,i) = kmeans(Vs,i,'emptyaction','singleton','replicate',5);
        end
        va=evalclusters(Vs,clust,'CalinskiHarabasz');
        optimk=va.OptimalK;
        
end
