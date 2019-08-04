
data=X';
thr=15;
L=25;
censoredornot=0; %=1 if you want censoring otherwise set it to 0
t=1;
msteps=10;
stepsize=0.2;


begin=0;
[logsigma,dim_norm] = diffusion_map_kt(data,censoredornot,thr,thr,thr+L,begin,msteps,stepsize);
[dim_max,id]=max(dim_norm); %find local maxima

sigma=10^(logsigma(id));
%%%%%%%%%%%%%%%%%%%%%%perform diffusion maps analysis%%%%%%%%%%%%%%%%%%
[psi,E] = diffusion_map_main(data,censoredornot,thr,thr,thr+L, t, sigma);
scatter3(psi(:,1),psi(:,2),psi(:,3),30);

