%%%% example in HiC networks
clc
clear
close all
addpath(genpath(pwd))
%%%load the example data. We only show the resolution 1K here.
% If you want to play with resolution 5K, you just need to change 1000 to
% 5000 when you do the loading as below

HiC_raw_resolution1K = load('mat_16_1000.txt'); % raw HiC contact network with resolution 1K
label_resolution1K = load('class_16_1000.txt'); % labels for each node in HiC network with resolution 1K

%%%enhance the network
HiC_NE=Network_Enhancement(HiC_raw_resolution1K);

% %%% run Louvain method to perform module detection/clustering

% IMPORTANT!!! please mex the two .cpp files first before you run Louvain method
% please run the following command in matlab: 
% mex jl_mnew.cpp
% mex jl_clust.cpp

com_raw = louvain(HiC_raw_resolution1K);
com_NE = louvain(HiC_NE);
% 
% %%% Calculate NMI between obtained clustering and ground-truth
NMI_raw = nmi(com_raw,label_resolution1K);
NMI_NE = nmi(com_NE,label_resolution1K);

fprintf('The NMI on raw HiC network is %6.4f \n', NMI_raw);
fprintf('The NMI on enhanced network is %6.4f \n', NMI_NE);


%% For visualization purposes, we first normalize the contact map and trim data trim data in the extreme 2% at each tail. 
% Values at the lower tail are set to zero and values at the higher tail are set to the 98% percentile



% Normalize 
% Normalize both matrix. The normalization factor is the square root of the
% column sum times row sum (see square root VC normalization in Rao 2014)


positions = [];
% fine the position where the label changes (boundary of domains)
for i =  1:18;
    positions = [positions, find(label_resolution1K==i,1)];
end
HiC_raw_resolution1K = log(1+HiC_raw_resolution1K);
HiC_NE = log(1+HiC_NE);
nHiC = sum(HiC_raw_resolution1K);
nNE = sum(HiC_NE);

for i = 1:size(HiC_raw_resolution1K,1) 
    for j = 1:size(HiC_raw_resolution1K,1)
        HiC_raw_resolution1K(i,j) = HiC_raw_resolution1K(i,j)/sqrt(nHiC(i)*nHiC(j));
        HiC_NE(i,j) = HiC_NE(i,j)/sqrt(nNE(i)*nNE(j));
    end
end

qminHiC = quantile(HiC_raw_resolution1K(HiC_raw_resolution1K>0), .02);
qminNE = quantile(HiC_NE(HiC_NE>0), .02);
qmaxHiC = quantile(HiC_raw_resolution1K(HiC_raw_resolution1K>0), .98);
qmaxNE = quantile(HiC_NE(HiC_NE>0), 0.98);
% 
HiC_raw_resolution1K(HiC_raw_resolution1K < qminHiC) = 0;
HiC_NE(HiC_NE<qminNE) = 0;
% 
HiC_raw_resolution1K(HiC_raw_resolution1K > qmaxHiC) = qmaxHiC;
HiC_NE(HiC_NE > qmaxNE) = qmaxNE;

% Plot figure 2A 
% Plots the first 10 domains from the NE treated data and the raw data. The
% lines correspond to where the boundary of the domains are (as they were
% stitched together)


ind = 1:positions(10);
clf
cmap = 'parula';
figure(1)
hold on
axes('Units', 'normalized', 'Position', [0 0 1 1])
imagesc(HiC_NE(ind,ind))
%cmap = cbrewer('seq','YlGn',100);
colormap(cmap)

for i = 1:10;
    line([positions(i+1) positions(i+1)], [positions(i) positions(i+1)], 'Color', 'r', 'LineWidth',4);
    line([positions(i) positions(i+1)], [positions(i) positions(i)], 'Color', 'r','LineWidth',4);
end
%line([1 500], [50 50], 'Color','r', 'linewidth', 1)
%plot(x,y, 'y', 'linewidth',5)
hold off
shg

figure(2)
hold on
axes('Units', 'normalized', 'Position', [0 0 1 1])
imagesc(HiC_raw_resolution1K(ind, ind))
for i = 1:10;
    line([positions(i+1) positions(i+1)], [positions(i) positions(i+1)], 'Color', 'r', 'LineWidth',4);
    line([positions(i) positions(i+1)], [positions(i) positions(i)], 'Color', 'r','LineWidth',4);
end
colormap(cmap)

