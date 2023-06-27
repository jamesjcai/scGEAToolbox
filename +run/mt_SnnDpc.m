function [c]=mt_SnnDpc(s,cluK,knnK)
%Clustering cell embeddings using SNNDPC - a SNN clustering algorithm
%
% http://mlwiki.org/index.php/SNN_Clustering#SSN_Clustering_Algorithm
% https://link.springer.com/article/10.1007/s12539-019-00357-4

if nargin<3, knnK=4; end
if nargin<2, cluK=10; end

pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'external','mt_SNNDPC');
if ~(ismcc || isdeployed)
    addpath(pth);
end
x=SnnDpc_ori(s,ones(size(s,1),1),knnK,'AutoPick',cluK);
c=x.cluster;

end