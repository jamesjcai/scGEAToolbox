function [s,c]=sc_viewngrpcells(Xcell,s)

X=[];
c=[];
for k=1:length(Xcell)
    X=[X Xcell{k}];
    n=size(Xcell{k},2);
    c=[c;k*ones(n,1)];
end
if nargin<2
    s=sc_tsne(X,3);
end
hold on
for k=1:length(Xcell)
    i=c==k;
    scatter3(s(i,1),s(i,2),s(i,3),10);
end

