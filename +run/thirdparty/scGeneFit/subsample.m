%Filename:      subsample.m
%Last edit:     Feb 11 2019
%Description:   samples a fraction of the data making sure all labels
%               are represented
%Inputs:
%               -data: 
%     
%               a G x C array of data points in dimension G.
% 
%               -labels: 
%
%               a 1 x C array labels between 1 and N
%
%               -fraction: 
%     
%               fraction of samples to be selected
%
% 
%Outputs:
%               -subsamples: 
% 
%               a G x c array with a random fraction of data points
%
%               -sublabels: 
% 
%               labels corresponding to the subsamples 
% 
%Documentation:
% 
% -------------------------------------------------------------------------

function [subsamples, sublabels]= subsample(data, labels, fraction)

[d,~]=size(data);
L=size(unique(labels),1);
count=zeros(L,1);
s=0;
for i=1:L
    count(i)=sum(labels==i);
    s=s+min(ceil(fraction* count(i))+2, count(i));
end

subsamples=zeros(d, s);
sublabels=zeros(s, 1);
idx=0;

for i=1:L
    aux=data(:,labels==i);
    auxl=labels(labels==i);
    N= min(count(i), ceil(fraction*count(i))+2); %at least 3 points per cluster
    I=randperm(size(aux,2), N);
    subsamples(:,idx+1:idx+N)=aux(:,I);
    sublabels(idx+1:idx+N)=auxl(I);
    idx=idx+N;
end

end
