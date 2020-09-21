function [ Matrix ] = MSNN( X,m,k )
% SNN calculates shared nearest neighbours based dissimilarity
% The SNN similarity of any two points is based on the overlap between
% their k-nearest neighbours lists. 

%   Input  : X  : data set (m,n); m-objects, n-variables 
%            k :  to build k-nearest neighbours lists
%            m :  original dissimilarity matrix (m by m)
%   Output : Matrix: SNN dissimilarity matrix (m by m)

I=zeros(k+1,size(m,1));
parfor i=1:size(m,1) 
    [~,index] = sort(m(i,:),2,'ascend');
    I(:,i)=index(1:k+1)';
end
     
    



n=size(X,1);
Matrix=zeros(n,n);
tic
m=triu(true(size(X,1)),1);
[posc,posr]=find(m==1);
parfor i=1:size(posc,1)
    if (size(find(I(2:end,posc(i))==posr(i)),1)+size(find(I(2:end,posr(i))==posc(i)),1))==2
        p(i)=size(intersect(I(2:end,posc(i)),I(2:end,posr(i))),1);
    else
        p(i)=0;
    end
end
for i=1:size(posc,1)
    Matrix(posc(i),posr(i))=p(i);
end
Matrix=Matrix'+Matrix;

Matrix=Matrix+eye(n,n)*k;
toc
disp('Matrix complete...');


end

