function [Eigen_Vector,Eigen_Value]=Find_K_Min_Eigen(Matrix, Eigen_NUM)  
  
[NN,NN]=size(Matrix);  
[V,S]=eig(Matrix); %Note this is equivalent to; [V,S]=eig(St,SL); also equivalent to [V,S]=eig(Sn,St); %  
  
S=diag(S);  
[S,index]=sort(S);  
  
Eigen_Vector=zeros(NN,Eigen_NUM);  
Eigen_Value=zeros(1,Eigen_NUM);  
  
p=1;  
for t=1:Eigen_NUM  
    Eigen_Vector(:,t)=V(:,index(p));  
    Eigen_Value(t)=S(p);  
    p=p+1;  
end  