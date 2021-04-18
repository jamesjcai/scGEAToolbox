function [W]=fun_sim_matrix(X,varargin)
% compute cell-to-cell similarity matrix
p = inputParser;
defaultType = 'simlr';
validTypes = {'simlr','soptsc','magic'};
checkType = @(x) any(validatestring(x,validTypes));

addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,X,varargin{:})

% pw1=fileparts(which(mfilename));
pw1 = fileparts(mfilename('fullpath'));
switch p.Results.type
    case 'simlr'        
        pth=fullfile(pw1,'thirdparty/SIMLR');
        addpath(pth);
        pth=fullfile(pw1,'thirdparty/SIMLR/src');
        addpath(pth);
        
    case 'magic'
        pth=fullfile(pw1,'thirdparty/MAGIC');
        addpath(pth);
        npca = 100;
        
        [X]=sc_norm(X,'type','libsize');
        X=log(X + 0.1);
        
        data=X';
        [U,~,~] = randPCA(data', npca); % this is svd
        pc = data * U; % this is PCA without mean centering to be able to handle sparse data
        % compute kernel
        disp 'computing kernel'
        W = compute_kernel(pc, 'k', 10, 'a', 15, 'distfun', 'euclidean');
        W=full(W);

        
            
    case 'soptsc'        
        pth=fullfile(pw1,'thirdparty/SoptSC');
        addpath(pth);
        pth=fullfile(pw1,'thirdparty/SoptSC/NNDSVD');
        addpath(pth);
        pth=fullfile(pw1,'thirdparty/SoptSC/symnmf2');
        addpath(pth);
        
        [X]=sc_norm(X,'type','libsize');
        X=log(X + 1);
                
        realdata = X;
        realdata = realdata-min(realdata(:));
        realdata = realdata./max(realdata(:));

        [~,n] = size(realdata);
        for i = 1:n
            realdata(:,i) = realdata(:,i)/norm(realdata(:,i));
        end
        lambda = 0.5;
        W=SimilarityM(realdata,lambda,X);        
end
