function W=censoring(data,thr,M1,M2,sigma) 
%thr = censoring value (e.g. 15)
%[M1,M2] = uncertainity range pos(e.g. 40-15)
n=size(data,1);
W=zeros(n,n);
KT=sigma^2;

M2=M2+sqrt(KT);
M1=M1-sqrt(KT);


for i=1:n
for j=1:n
 
          W(i,j)=1;
         for g=1:size(data,2)
               if (data(i,g)==thr)
                  
               if data(j,g)==thr
                  W(i,j)=W(i,j);
               else
                  W(i,j)= W(i,j)*(pi*KT/2)^(-1/4)*(pi*KT/4)^(1/2)*( erfc( (M1-data(j,g))/sqrt(KT)) - erfc( (M2-data(j,g))/sqrt(KT)) )/sqrt(M2-M1);
                 
               end
               elseif (data(j,g)==thr) 
                  W(i,j)= W(i,j)*(pi*KT/2)^(-1/4)*(pi*KT/4)^(1/2)*( erfc( (M1-data(i,g))/sqrt(KT)) - erfc( (M2-data(i,g))/sqrt(KT)) )/sqrt(M2-M1);

               else
         
                  W(i,j)= W(i,j)* exp(-(data(i,g)-data(j,g))^2/(KT*2));
              end
              
         end
end
end


