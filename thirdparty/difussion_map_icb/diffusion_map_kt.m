function [logsigma,dim_norm] = diffusion_map_kt(data,censored,thr,M1,M2,begin,m,stepsize)

[n,G]=size(data);
avrdnorm=zeros(m,1);
logsigma= zeros(m,1);

for Dround=1:m
      sigma=10^(begin+Dround*(stepsize));

if censored==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%censoring%%%%%%%%%%%%%%%%%%%%%%%%%
W=censoring(data,thr,M1,M2,sigma); 
else
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%no censoring%%%%%%%%%%%%%%%%%%%%%%
d2=zeros(n,n);
for g=1:G 
    d2_g = (bsxfun(@minus,reshape(data(:,g),n,1),reshape(data(:,g),1,n))).^2;  
    d2=d2+d2_g;
end
W=exp(-d2/(2*sigma^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

diagD=nansum(W);

avrdnorm(Dround)=(sum(log10(diagD/n)./(diagD)))/sum(1./diagD);
logsigma(Dround)= log10(sigma);
end
 
dim_norm=diff((avrdnorm))./diff(logsigma);

