function [c]=sc_louvain(A)
%Louvain clustering algorithm with adjacency matrix A

pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/Network_Enhancement/examples/HiC_Network');
addpath(pth);
c = louvain(A);


