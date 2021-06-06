function [a,fano,cv]=i_grpmean(x,c)

n=numel(unique(c));
a=zeros(size(x,1),n);
%atrim=zeros(size(x,1),n);
b=zeros(size(x,1),n);
c=zeros(size(x,1),n);
for k=1:n
%    atrim(:,k)=full(trimmean(x(:,c==k),10,2));
    a(:,k)=full(mean(x(:,c==k),2));
    b(:,k)=full(var(x(:,c==k),0,2));
    c(:,k)=full(std(x(:,c==k),0,2));
end
if nargout>1
    fano=b./a;
end
if nargout>2
    cv=c./a;
end
%if nargout==1
%    a=atrim;
%end
end
