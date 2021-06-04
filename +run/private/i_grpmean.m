function a=i_grpmean(x,c)
n=numel(unique(c));
a=zeros(size(x,1),n);
for k=1:n
    a(:,k)=full(mean(x(:,c==k),2));
end
end