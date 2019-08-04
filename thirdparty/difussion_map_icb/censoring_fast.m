function [W]=censoring_fast(data,thr,M1,M2,sigma)

M=max(thr,thr+L);
thr=min(thr,thr+L);
KT=sigma^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exp(-a*(X-XA)^2)*exp(-b*(X-XB)^2)
a=(1/KT);
b=(1/(sqrt(M2-M1)*KT));
%XA is the center of the normal Gaussian, hence data(i,j) 
XB=(M1+M2)/2;   %XB is the center of the censoring Gaussian

norm_ratio=sqrt((a+b)/(2*a))*(b/a)^(1/4); % normalization coefficients ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=size(data,1);
W=zeros(n,n);
count=zeros(n,n);

for i=1:n
for j=1:n
          count(i,j)=0;  
          W(i,j)=1;   
         for g=1:size(data,2)
               if (data(i,g)==thr)
               if data(j,g)==thr
                  W(i,j)=W(i,j)*sqrt(1);
               else
                    W(i,j)= W(i,j)+(a*b/(a+b))*(data(j,g)-XB)^2 + log (norm_ratio);
                 count(i,j)=count(i,j)+1;
               end
               elseif (data(j,g)==thr) 
                  W(i,j)= W(i,j)+ (a*b/(a+b))*(data(i,g)-XB)^2 + log (norm_ratio);
               else
%         
 W(i,j)= W(i,j)+ ((a/2)*(data(i,g)-data(j,g))^2);%
              end
              
         end
end
end
W=exp(-W);

