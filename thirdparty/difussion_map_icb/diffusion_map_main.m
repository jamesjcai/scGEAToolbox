 function [psi,E,D] = diffusion_map_main(data,censored,thr,M1,M2, t, sigma)

%data=high dimensional matrix data to be analysied
%set censored=1 if censoring required, if not set censored=0
%thr =  censoring values fro each gene(e.g. 15)
%[M1,M2] = vecotr of uncertainity range pos(e.g. [15,40])
%t= 1 for with density normalisation factor, t=0 for without density
%normalisation
%sigma=diffusion scale parameter of the Gauusian kernel
% psi= mapped coordinates of each cell
% E= first (largest) eigenvalues of the transition matrix, a gap in E shows indicates the suitable no. of mapping dims
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[n,G]=size(data);
if censored==1;
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

D=sum(W);
q=(D'*D).^t;
W(1:n+1:end)=0; %zero diagonal


H=W./q;  %density normalisation

D_=diag(sum(H));

Hp=D_^(-1)*H;
 
[psi_nsort,En]=eig(Hp);
%[psi_nsort,En]=eigs(Hp); %to calculate only few largest
%eigenvalues/eigenvectors


[E, indE]=sort(real(diag(En)),'descend'); 
psi=psi_nsort(:,indE); 

psi(:,1)=[];
E(1)=[];

toc/60
    
