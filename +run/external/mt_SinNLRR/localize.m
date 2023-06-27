function [localX,coverage] = localize( C )
%C is the coefficient matrix
%[tmp,ind] = sort(C,1,'descend');
[m,n]=size(C);
localX=C;
coverage=zeros(1,n);
for i=1:n
   thr=C(i,i)/1.5;
   localX(localX(:,i)<thr,i)=0;
   coverage(1,i)=mean(C(i,i)./localX(localX(:,i)>thr,i));
end
end

