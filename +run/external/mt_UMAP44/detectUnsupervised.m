function [cnt, unsupervisedIdxs]=detectUnsupervised(umap, inData, ...
    sduLimit, parameterLimit)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

supervisors=umap.supervisors;
subsetIds=unique(supervisors.labels);
subsetIds=subsetIds(subsetIds>0);
nSubsets=length(subsetIds);

a=zeros(nSubsets, size(inData,1));
for i=1:nSubsets
    r=umap.raw_data(supervisors.labels==subsetIds(i),:);
    means_=mean(r);
    stds_=std(r);
    B=(abs(inData-means_(1,:)))./stds_(1,:);
    a(i,:)=(sum(B>sduLimit, 2)>parameterLimit)';
end
S=sum(a);
unsupervisedIdxs=S==nSubsets;
cnt=sum(unsupervisedIdxs);
end
        