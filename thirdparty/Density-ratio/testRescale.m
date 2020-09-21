%% test ReScale
% this script compares the performance of original clustering algorithms
% with the performance of their ReScale versions. Parameters for each
% algorithm are tuned to the best. The results are the same as Table 6 of
% the original paper.

clear
clc
load('S1');  % load a dataset

%% Original version
%%%%%%%%% DBSCAN
 
SimMatrix=pdist2(data,data,'minkowski',2);
eps=0.02;
MinPts=3;
Tclass=Mdbscan(data,MinPts,eps,SimMatrix)';  
[ ~,~,~,~,~,~,~,~,DBSCAN_Fmeasure,~] = evaluate(class,Tclass) % the best performance

%%%%%%%%% SNN
k=ceil(size(data,1)^0.5);
SimMatrix= SNN(data,k);
SimMatrix=k-SimMatrix;
Eps=10;
MinPts=9;
Tclass=Mdbscan(data,MinPts,Eps,SimMatrix)'; 
[ ~,~,~,~,~,~,~,~,SNN_Fmeasure,~] = evaluate(class,Tclass) % the best performance 

%% ReScale version
%%%%%%%%%%% DBSCAN
psi=10;
eta=0.1;
[ ndata ] = Rescale( psi,eta,data);
 
SimMatrix=pdist2(ndata,ndata,'minkowski',2);
eps=0.218;
MinPts=10;
Tclass=Mdbscan(ndata,MinPts,eps,SimMatrix)'; 
[ ~,~,~,~,~,~,~,~,ReScale_DBSCAN_Fmeasure,~] = evaluate(class,Tclass) % the best performance

%%%%%%%%%%%%%% SNN
psi=1000;
eta=0.5;
[ ndata ] = Rescale( psi,eta,data);

k=ceil(size(ndata,1)^0.5);
SimMatrix= SNN(ndata,k);
SimMatrix=k-SimMatrix;
Eps=11;
MinPts=10;
Tclass=Mdbscan(ndata,MinPts,Eps,SimMatrix)'; % SNN
[ ~,~,~,~,~,~,~,~,ReScale_SNN_Fmeasure,~] = evaluate(class,Tclass) % the best performance
 