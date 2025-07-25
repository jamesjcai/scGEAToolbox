function [R] = get_R(Z,Y,sigma)
%Z2=L*Z;
N=size(Z,2);
K=size(Y,2);
R=[];
for i=1:N
    sum=0;
    for k=1:K
        sum=sum+exp(-((Z(:,i)-Y(:,k))'*(Z(:,i)-Y(:,k)))/sigma);
    end
    
    for j=1:K
        R(i,j)=exp(-((Z(:,i)-Y(:,j))'*(Z(:,i)-Y(:,j)))/sigma)/sum;
    end
end
end

