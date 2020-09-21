% This script compares the MDS plot of a original data with its rescaled
% version. MDS is a data analysis technique for visualising the information contained 
% in a dissimilarity matrix. The aim of an MDS algorithm is to place each data point 
% in a p-dimensional space (where p is the desired number of dimensions, specified by 
% a user), while preserving as well as possible the pairwise dissimilarities between 
% objects. Note that the orientation of the reconstruction is arbitrary when using 
% MDS since the transformation aims to preserve the pairwise dissimilarities.
% Figure 10 of the original paper is based on this script.

load('hard') % load a dataset

%% Visualise original data
figure
if size(data,2)==2
    gscatter(data(:,1),data(:,2),class,'rgbmykc','xdo+ps*')
    box off
    set(gcf,'color','w');
else
    Matrix=pdist2(data,data,'minkowski',2);    %% plot MDS if ndata has more than 2 attributes
    Y = mdscale(Matrix,2,'criterion','stress','start','random');
    gscatter(Y(:,1),Y(:,2),class,'rgbmykc','xdo+ps*')
    box off
    set(gcf,'color','w');
end
title('Original data')
%% ReScale
psi=100
eta=0.1
[ ndata ] = Rescale( psi,eta,data);


%% Visualise scaled data
figure
if size(data,2)==2
    gscatter(ndata(:,1),ndata(:,2),class,'rgbmykc','xdo+ps*')
    box off
    set(gcf,'color','w');
else
    Matrix=pdist2(ndata,ndata,'minkowski',2);    %% plot MDS if ndata has more than 2 attributes
    Y = mdscale(Matrix,2,'criterion','stress','start','random');
    gscatter(Y(:,1),Y(:,2),class,'rgbmykc','xdo+ps*')
    box off
    set(gcf,'color','w');
end
title('ReScaled data')