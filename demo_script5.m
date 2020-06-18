%% Demonstration of Pseudotime Analysis and Gene Network Functions in scGEAToolbox
%% Load example data set, X
%
cdgea; % set working directory
[X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');

%% Select genes with at least 3 cells having more than 5 reads per cell. 
%
[X,genelist]=sc_selectg(X,genelist,5,3);

%% Trajectory analysis using the PHATE+splinefit method
%
% s=run_phate(X,3,false,true);
% [t,xyz1]=i_pseudotime_by_splinefit(s,1);
% hold on
% plot3(xyz1(:,1),xyz1(:,2),xyz1(:,3),'-r','linewidth',2);

% Calculte pseudotime T
figure;
t=sc_trajectory(X,"type","splinefit","plotit",true);

%% Plot gene expression profile of cells ordered according to their pseudotime T.
%

r=corr(t,X','type','spearman'); % Calculate linear correlation between gene expression profile and T
[~,idxp]= maxk(r,4);  % Select top 4 positively correlated genes
[~,idxn]= mink(r,3);  % Select top 3 negatively correlated genes
selectedg=genelist([idxp idxn]);

% Plot expression profile of the 5 selected genes
figure;
i_plot_pseudotimeseries(log(1+X),genelist,t,selectedg)

% % Nonlinear correlation  
%
% r=zeros(size(X,1),1);
% for k=1:size(X,1)
%     k
%     r(k)=distcorr(t,double(X(k,:))');
% end
% [~,idxp]= maxk(r,3);
% [~,idxn]= mink(r,2);
% selectedg=genelist([idxp; idxn]);
% figure;
% i_plot_pseudotimeseries(log(1+X),genelist,t,selectedg)

%% Trajectory analysis using TSCAN
%
% Calculte pseudotime T
figure;
t=sc_trajectory(X,"type","tscan","plotit",true);

r=corr(t,X','type','spearman'); % Calculate linear correlation between gene expression profile and T
[~,idxp]= maxk(r,4);  % Select top 4 positively correlated genes
[~,idxn]= mink(r,3);  % Select top 3 negatively correlated genes
selectedg=genelist([idxp idxn]);

% Plot expression profile of the 5 selected genes
figure;
i_plot_pseudotimeseries(log(1+X),genelist,t,selectedg)

%% Construct single-cell gene regulatory network (scGRN)

%% Using principal component regression (PCNet) method
%
X50=X(1:50,:);
genelist50=genelist(1:50);
A=sc_pcnet(X50);

% Plot constructed network
%
A=A.*(abs(A)>quantile(abs(A(:)),0.9));
G=digraph(A,genelist50);
LWidths=abs(5*G.Edges.Weight/max(G.Edges.Weight));
LWidths(LWidths==0)=1e-5;
figure;
plot(G,'LineWidth',LWidths);
p.MarkerSize = 7;
p.Marker = 's';
p.NodeColor = 'r';

%% Using GENIE3 method
%
X20=X(1:20,:);
genelist20=genelist(1:20);
A=run_genie3(X20);

% Plot constructed network
%
A=A.*(abs(A)>quantile(abs(A(:)),0.9));
G=digraph(A,genelist20);
LWidths=abs(5*G.Edges.Weight/max(G.Edges.Weight));
LWidths(LWidths==0)=1e-5;
figure;
plot(G,'LineWidth',LWidths);
p.MarkerSize = 7;
p.Marker = 's';
p.NodeColor = 'r';

%% The End