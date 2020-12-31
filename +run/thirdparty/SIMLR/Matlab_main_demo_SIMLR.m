% clear
% clc
% close all
% 
% addpath('data')
% addpath('src')
% dataset = {'1_mECS', '2_Kolod', '3_Pollen', '4_Usoskin'}
% 
% for i = 1:4
%     
%     % perform the analysis for the current dataset
%     load(['Test_' dataset{i}]);
%     C = max(true_labs); %%% number of clusters
%     rng(i,'twister'); %%% for reproducibility

load ../../example_data/example10xdata.mat

in_X=X';
C=5;
    [y, S, F, ydata,alpha] = SIMLR(in_X,C,10);
    
%     % report NMI values
%     NMI_i = Cal_NMI(y,true_labs);
%     fprintf(['The NMI value for dataset ' dataset{i} ' is %f\n'], NMI_i);
    
    % visualization
    figure;
   % gscatter(ydata(:,1),ydata(:,2),true_labs);
    gscatter(ydata(:,1),ydata(:,2),y);
    
%     ydata=sc_tsne(X,2,false,false);  % s=sc_tsne(X,ndim,plotit,donorm,dolog1p
%     gscatter(ydata(:,1),ydata(:,2),y);
% end
