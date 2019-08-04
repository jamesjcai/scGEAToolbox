function [S,c,cindex,D,DX]=distortionMinimization(X,phi0,k,c_0,DspFlag,epsilon);

% function [S,c,cindex,D,DX]=distortionMinimization(X,phi0,k,c_0,DspFlag,epsilon);
%
% Input:
%   X(n,d)        --- diffusion coordinates (right eigenvecs
%                     rescaled with eigenvals); each row is a data pt
%   phi0(n,1)     --- stationary distribution (first left
%                     eigenvector of Markov matrix) 
%   k             --- number of clusters
%   c_0(k,d)      --- initial centers
%   DspFlag       --- display flag
%   epsilon       --- relative distortion (stopping criteria)
%
% Output:
%   S(n,1)        --- labeling; n-dim vector with numbers between 1 and k
%   c(k,d)        --- geometric centroids
%   cindex(k,1)   --- diffusion centers -- subset of original data; 
%                     k-dim vector (row indices of X)
%   D             --- distortion
%   DX(n,k)       --- squared distance from each point to every centroid
%
% Original version by Stephane Lafon, Yale University, April 2005.
% Last modified, Ann B. Lee, 4/11/05. ABL at CMU: 5/7/08 

global cPoints  % Warning: global variable

[n,d]=size(X);
if(nargin<6)
    epsilon=1e-3;
end
if(nargin<5)
    DspFlag=0;
end
if(nargin<4)
    tmp_ind = ceil(rand(k,1)*n);  % random subset of X
    c_0 = X(tmp_ind,:);   % k-by-d matrix
    %S_0=ceil(rand(n,1)*k);
end
col='rgkbmcy';
%S=S_0;
%c=zeros(k,d);
c=c_0;
oldD=Inf;

MaxIter=1000;

for i=1:MaxIter,  % KMEANS LOOP
        
    %-----------------------------------
    % K-MEANS
    %-----------------------------------
    % Update distances to centroids and labels of data points:
    DX=[];    
    for j=1:k,
        dX=X-ones(n,1)*c(j,:);  % n-by-d
        DX=[DX dot(dX',dX')'];  % enter column wise; n-by-k matrix of distances
    end
    [Dtmp,j]=min(DX');  % min(DX,2)  -> 1-by-n
    S=j'; % new labels
    %D=Dtmp*phi0; % distortion (a number)
    % Check for empty clusters:
    for j=1:k,
        ind=find(S==j);
        if(isempty(ind)) % if cluster j is empty,
            [mx,m]=max(Dtmp); % find data point that is furthest from its centroid
            S(m)=j; % assign this point the label j
            Dtmp(m)=0;
            %c(j,:) = X(ind(m),:);  % make this point centroid j (redundant)
        end
    end    
    
    % Update centroids:
    for j=1:k,
        ind=find(S==j);
        c(j,:)=phi0(ind,1)'*X(ind,:)/sum(phi0(ind,1)); % find centroid of cluster j
        %plot(c(j,1),c(j,2),'k*');
    end
    
    % Distortion
    D=Dtmp*phi0;
    
    % Plot results:
    if DspFlag
        figure(1), clf, scatter(cPoints(:,1),cPoints(:,2),15,S); axis image, hold on;
        for j=1:k,
            dX=X-ones(n,1)*c(j,:);  % n-by-d
            DX=dot(dX',dX')';  % n-dim vector
            [dummy,tmpind]=min(DX,[],1);
            cindex(j)=tmpind;
            plot(cPoints(cindex(j),1),cPoints(cindex(j),2),'rx'); % diffusion centers
        end
        hold off;
        pause,
    end
     
    % Stopping criteria:
    if((oldD-D)/D < epsilon)
        break;
    end
    oldD=D;
end

%------------------------------------------------------------
% centroids => "diffusion centers" (subset of original data)
%------------------------------------------------------------
cindex=zeros(k,1);
for j=1:k,
    %ind=find(S==j);
    %dX=X(ind,:)-ones(length(ind),1)*c(j,:);
    dX=X-ones(n,1)*c(j,:);  % n-by-d
    DX=dot(dX',dX')';  % n-dim vector
    [dummy,tmpind]=min(DX,[],1);
    cindex(j)=tmpind;
    %DX=dot(dX',dX')';
    %[tmp,i]=min(DX);
    %cindex(j)=ind(i);
end
