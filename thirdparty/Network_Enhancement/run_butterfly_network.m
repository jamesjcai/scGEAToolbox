clear
clc
close all
addpath(genpath(pwd));

%%% load the raw network

load('Raw_butterfly_network.mat')

%%%run Network_Enhancement
W_butterfly_NE=Network_Enhancement(W_butterfly0);


%% print/plot the results
[~,acc_raw] = CalACC(W_butterfly0, labels); % calculate acc on the raw network

[~,acc_NE] = CalACC(W_butterfly_NE, labels); % calculate acc on the denoised network

fprintf('The accuracy on raw network is %6.4f \n', acc_raw);
fprintf('The accuracy on enhanced network is %6.4f \n', acc_NE);


figure;
NUM = 80; %the number of images per class
[ tpr0 ] = cal_specific_accuracy(W_butterfly0,labels,NUM);
plot((1:NUM), (tpr0), 'b-', 'Linewidth',5,'MarkerSize',5); hold on;
[ tpr1 ] = cal_specific_accuracy(W_butterfly_NE,labels,NUM);
plot((1:NUM), (tpr1), 'r-', 'Linewidth',5,'MarkerSize',5); hold on;
axis([0,80,0.45,1])
legend('Raw', 'NE');

h = xlabel('Number of Retrieval');set(h,'FontSize',16);
h = ylabel('Identification Accuracy');set(h,'FontSize',16);

set(gca,'FontSize',16)







