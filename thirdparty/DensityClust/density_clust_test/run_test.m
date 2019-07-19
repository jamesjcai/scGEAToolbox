load s1131_cr.mat
X=t_sne;
%% Settings of System Parameters for DensityClust
dist = pdist2(X, X); % [NE, NE] matrix (this case may be not suitable for large-scale data sets)
% average percentage of neighbours, ranging from [0, 1]
% as a rule of thumb, set to around 1%-2% of NE (see the corresponding *Science* paper for more details)
percNeigh = 0.02;
% 'Gauss' denotes the use of Gauss Kernel to compute rho, and
% 'Cut-off' denotes the use of Cut-off Kernel.
% For large-scale data sets, 'Cut-off' is preferable owing to computational efficiency,
% otherwise, 'Gauss' is preferable in the case of small samples (especially with noises).
kernel = 'Gauss';
% set critical system parameters for DensityClust
[dc, rho] = paraSet(dist, percNeigh, kernel); 
figure(2);
plot(rho, 'b*');
xlabel('ALL Data Points');
ylabel('\rho');
title('Distribution Plot of \rho');

%% Density Clustering
isHalo = 1; 
[numClust, clustInd, centInd, haloInd] = densityClust(dist, dc, rho, isHalo);

figure;
i_myscatter(X,clustInd)
