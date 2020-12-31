function [c,cls]=sc_sc3(X,k,plotit)
% SC3 - consensus clustering of single-cell RNA-seq data
% Ref: https://www.nature.com/articles/nmeth.4236

if nargin<3
    plotit=false;
end
if nargin<2
    % [optimk]=fun_num_cluster(X,'type','simlr');
    % [optimk]=fun_num_cluster(X,'type','sc3');
    [optimk]=fun_num_clusters(X);
else
    optimk=k;
end

[X]=sc_norm(X);
X=log2(X+1);

drange=get_drange(X);

Dis=squareform(pdist(X'));
[cls1]=get_clusterarray(Dis,optimk,drange);

Dis=1-corr(X,'type','s');
[cls2]=get_clusterarray(Dis,optimk,drange);

Dis=1-corr(X,'type','p');
[cls3]=get_clusterarray(Dis,optimk,drange);

oldpath=pwd;
pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/ClusterPack');
addpath(pth);
cd(pth);
cls=[cls1;cls2;cls3];
c=clusterensemble(cls,optimk);
cd(oldpath);

if plotit
    clusion(Dis,c);
end

end



function [cls]=get_clusterarray(Dis,optimk,drange)
    [Vs1]=pca(Dis);
    [Vs2]=transform_Laplacian(Dis,max(drange));
    cls=[];
    for j=1:length(drange)
        idx=kmeans(Vs1(:,1:drange(j)),optimk,'MaxIter',1e9,'emptyaction','singleton','replicate',5);
        cls=[cls; idx'];
        idx=kmeans(Vs2(:,1:drange(j)),optimk,'MaxIter',1e9,'emptyaction','singleton','replicate',5);
        cls=[cls; idx'];
    end
end

function drange=get_drange(X)
    % https://www.nature.com/articles/nmeth.4236
    % the best clusterings were achieved when d was between 4% and 7% of the number of cells, N (Fig. 1c, Supplementary Fig. 3a and Online Methods).
    n=size(X,2);
    drge=round(n.*[0.04 0.07]);
    drange=drge(1):drge(2);
    if length(drange)>15
        dx=drange(randperm(length(drange)));
        drange=dx(1:15);    
    end
end

function [V]=transform_Laplacian(Dis,k)

    A=exp(-Dis./max(Dis(:)));   % adjacency matrix
%     xD=diag(sum(A).^-0.5);  % D=diag(sum(A)); % d(i) the degree of node i
%     xA=xD*A*xD;             % normalized adjacenty matrix
%     L=eye(size(A,1))-xA;    % also L=xD*(D-A)*xD 
    [~,L]=i_sbe_laplacian_matrix(A);

    % see https://people.orie.cornell.edu/dpw/orie6334/lecture7.pdf
    % see https://en.wikipedia.org/wiki/Laplacian_matrix#Symmetric_normalized_Laplacian_2

%     [V,D]=eig(L);
%     [~,ind]=sort(diag(D));
%     V = V(:,ind);
      [V,~]=eigs(L,k,'smallestreal');
end


function [L,Lnorm] = i_sbe_laplacian_matrix(A)
% Get graph Laplacian matrix
%
%   L = laplacian(g)
%
% graph Laplacian matrix is defined by L = D - A, where D is vertex degree
% diagonal matrix and A is adjacency matrix.
%
% See also: adjacency
% L = diag(sum(A)) - A;

% Systems Biology & Evolution Toolbox
% Author: James Cai
% Email: jcai@tamu.edu
% Website: https://github.com/jamesjcai/SBEToolbox_lite

% https://github.com/dtuia/KEMA/blob/7378c0fce50a818c2fb59f5de9344ea5c1929fa4/general_routine/laplacian.m
% https://github.com/KavehFathian/clear/blob/e72318ad52082485442f55ad5a7ee9b989b85677/Algorithms/Helpers/NormalizeLap.m
% https://github.com/hungrydoggy/Pinocchio/blob/5664503b210005fa4f1fc053e237ba3ecf6a7945/skeletonizer/matlab/toolbox/compute_mesh_laplacian.m
% https://github.com/fljohnston1/otto-group-product/blob/02a6f35f8c144ed52a2097a1965d561172f2e701/SpectralClustering.m

D=sum(A);
L=diag(D)-A;
if nargout>1
    D(D~=0)=sqrt(1./D(D~=0));
    D=diag(D);
    % Lnorm=D*L*D;
    Lnorm=eye(size(A,1))-D*A*D;     % L = I-D^-1/2*W*D^-1/2
end
end