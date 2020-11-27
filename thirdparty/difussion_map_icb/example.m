
dataStruc=load('data/Guo_data.mat');
data=dataStruc.data;
Features=load('data/Guo_features.mat');
labels=Features.labels; 
genes=Features.genes;
thr=Features.thr;
L=Features.L;

censoredornot=0; %=1 if you want censoring otherwise set it to 0
t=1;
msteps=10;
stepsize=0.2;

%%
%%%%%%%%%%%%%%%this part is to estimate a proper sigma, not needed if a suitable sigma is already deceided%%%
%with censoring this part can take about 7 mins on typical PCs. Without
%censoring much faster (about one second)
begin=0;
[logsigma,dim_norm] = diffusion_map_kt(data,censoredornot,thr,thr,thr+L,begin,msteps,stepsize);
[dim_max,id]=max(dim_norm); %find local maxima

figure;
plot(logsigma(1:(msteps-1)),dim_norm(1:(msteps-1)),'g*-','LineWidth',1);

sigma=10^(logsigma(id));
%%%%%%%%%%%%%%%%%%%%%%perform diffusion maps analysis%%%%%%%%%%%%%%%%%%
[psi,E] = diffusion_map_main(data,censoredornot,thr,thr,thr+L, t, sigma);
%%%%%%%%%%%%%%%plotting eigenvalues and the map on the largest eigenvectors (%%%%%%%%%%%%%%%%%%%
figure;
plot(E(1:19),'k*-','LineWidth',1);
xlabel('i','FontSize',20);
ylabel('\lambda_{i}','FontSize',20);
figure;
%randomize cells order for plotting such that one colour (corresponding to a label) is not always on top of one other%%%%
n=size(data,1);
ix = randperm(n);
Ashuffled = psi(ix,:);
labTshuffled=labels(ix);

scatter3(Ashuffled(:,1),Ashuffled(:,2),Ashuffled(:,3),30,labTshuffled,'fill');

set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'ZTick',[])
xlabel('DC1','FontSize',20);
ylabel('DC2','FontSize',20);
zlabel('DC3','FontSize',20);