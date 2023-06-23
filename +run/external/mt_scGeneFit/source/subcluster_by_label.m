%Filename:      subcluster_by_label.m
%Last edit:     Feb 11 2019
%Description:   selects the indices of a subcluster
%               
%Inputs:
%               -labels1: 
%     
%               a 1 x C array of labels.
% 
%               -val_label: 
%
%               int indicating the labels1 to select
%
%               -labels2: 
%     
%               labels of the next level of the hierarchy
%
% 
%Outputs:
%               -indices: 
% 
%               indices of the data points where labels1==val_label 
%
%               -sublabels: 
% 
%               labels corresponding to next level hierarchy 
% 
%Documentation:
% 
% -------------------------------------------------------------------------
function [indices,labels]= subcluster_by_label(labels1, val_label, labels2)

indices=labels1==val_label;

auxlabels=labels2(indices);
values=unique(auxlabels);
labels=zeros(size(auxlabels));

for i=1:size(values,1)
    labels(auxlabels==values(i))=i;
end

end
