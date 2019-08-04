% Example code for diffusion maps (diffuse.m) 
% and diffusion coarse-graining (diffusion_kmeans.m)
%
% Ann Lee, March 2009. Last modified, JWR: 3/23/2009

%--------------------------------------------------------------------
% LOAD DATA
%--------------------------------------------------------------------
close('all');  % close all figures
load('annulus.mat'); 
[n,p]=size(Data); % Data(n,p), where n=#observations, p=#variables
plot(Data(:,1),Data(:,2),'.'); axis image
title('Example data');
D=squareform(pdist(Data)); % pairwise distances, n-by-n matrix

%--------------------------------------------------------------------
% CHECK THE DISTRIBUTION OF NEAREST NEIGHBORS
%--------------------------------------------------------------------
test_flag=1;
if test_flag % Estimate distribution of k-nearest neighbor
    D_sort = sort(D,2);
    k=30, %30;
    dist_knn = D_sort(:,1+k);  % distance to k-nn
    median_val = median(dist_knn), eps_val = median_val^2/2,
    sigmaK = sqrt(2)*median_val;
    figure, hist(dist_knn); colormap cool;
    title('Distribution of distance to k-NN');
end

%--------------------------------------------------------------------
% SET PARAMETERS IN MODEL
%--------------------------------------------------------------------
eps_val=0.05;  
neigen=10;
flag_t=0; %flag_t=0 => Default: multi-scale geometry
if flag_t
    t=3;  % fixed time scale  
end

%--------------------------------------------------------------------
% EIGENDECOMPOSITION
%--------------------------------------------------------------------
% Call function "diffuse.m"
[X, eigenvals, psi, phi] = diffuse(D,eps_val,neigen);


%--------------------------------------------------------------------
% DIFFUSION K-MEANS
%--------------------------------------------------------------------
k=2;  % number of clusters
Niter=100;
phi0=phi(:,1);
[idx, C, ERR, DX] = diffusion_kmeans(X, k, phi0, Niter);

%--------------------------------------------------------------------
% PLOT RESULTS
%--------------------------------------------------------------------
figure, % Fall-off of eigenvalues
subplot(2,1,1), plot(eigenvals(2:neigen+1),'*-'); % non-trivial eigenvals
ylabel('Eigenvalues \lambda');
if flag_t
    lambda_t=eigenvals(2:neigen+1).^t; 
    subplot(2,1,2), plot(lambda_t,'*-');  % multiscale weighting
    ylabel('\lambda^t');  
else
    lambda_multi = eigenvals(2:neigen+1)./(1-eigenvals(2:neigen+1));
    subplot(2,1,2), plot(lambda_multi,'*-');  % multiscale weighting
    ylabel('\lambda / (1-\lambda)');
end
figure, % Diffusion map
scatter3(X(:,1),X(:,2),X(:,3),10,'b'); 
title('Embedding with first 3 diffusion coordinates');
xlabel('X_1'); ylabel('X_2'); zlabel('X_3');
hold on
scatter3(C(:,1),C(:,2),C(:,3),60,'k','filled'); 

figure,  % K-means labeling
scatter(Data(:,1),Data(:,2),20,idx,'filled'); axis image
title('K-means with K=2');
