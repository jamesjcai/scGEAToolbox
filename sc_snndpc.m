function [c]=sc_snndpc(s,cluK,knnK)
%Clustering cell embeddings using SNNDPC - a SNN clustering algorithm
%
% http://mlwiki.org/index.php/SNN_Clustering#SSN_Clustering_Algorithm
% https://link.springer.com/article/10.1007/s12539-019-00357-4
%
%
% https://doi.org/10.1016/j.ins.2018.03.031
% https://github.com/liurui39660/SNNDPC/tree/MatlabImplementation

if nargin<3, knnK=4; end
if nargin<2, cluK=10; end
c=run.SnnDpc(s,cluK,knnK);

end
