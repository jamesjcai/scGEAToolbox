function i_plot_pseudotimeseries(X,genelist,t,genes)
%Plot pseudotime series
if nargin<4
    genes=string(['AICDA','BACH2','BCL6','IRF4','PAX5','PRDM1','REL','RELA']);
end
[t,i]=sort(t);
n=length(genes);

% idx=zeros(n,1);
% for k=1:n    
%     ixx=find(genelist==genes(k));
%     if ~isempty(ixx)
%         idx(k)=find(genelist==genes(k));
%     end
% end

idx=[];
glist=[];
for k=1:n    
    ixx=find(genelist==genes(k));
    if ~isempty(ixx)
        idx=[idx; ixx];
        glist=[glist; genes(k)];
    end
end
x=X(idx,i);



co=[0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    0   0   1];
co=[co;co;co;co];
set(groot,'defaultAxesColorOrder',co(1:length(idx),:))


% figure;
hold on
msz=1;
for k=1:size(x,1)
    plot(t,x(k,:),'.','markersize',msz); 
end
xlabel('Pseudotime')
ylabel('Expression')

pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'+run','thirdparty','locfit','m');
addpath(pth);
pth=fullfile(pw1,'+run','thirdparty','locfit','mex');
addpath(pth);
pth=fullfile(pw1,'+run','thirdparty','locfit','source');
addpath(pth);

if size(t,2)~=1
    t=t';
end
Pk=[];
for k=1:size(x,1)
    fitm1 = locfit(t,x(k,:)');
    Pk=[Pk plot(t,predict(fitm1,t),'-','LineWidth',3)];
end
box on
% legend([kglist glist],'location','eastoutside')
legend(Pk,glist,'location','eastoutside')
xlim([0 max(t)])

end