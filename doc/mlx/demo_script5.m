%[text] # Demo 5 - Pseudotime Analysis and Gene Network Functions
%[text] ## Load example data set, X
cdgea; % set working directory
[X,genelist]=sc_readfile(pkg.i_exampledata('GSM3044891_GeneExp.UMIs.10X1.txt'));
%%
%[text] ## Select genes with at least 3 cells having more than 5 reads per cell.
[X,genelist]=sc_selectg(X,genelist,5,3);
%%
%[text] ## Trajectory analysis using the PHATE+splinefit method
%[text] s=run\_phate(X,3,true); \[t,xyz1\]=i\_pseudotime\_by\_splinefit(s,1); hold on plot3(xyz1(:,1),xyz1(:,2),xyz1(:,3),'-r','linewidth',2);
% Calculte pseudotime T
figure;
t=sc_trajectory(X,"type","splinefit","plotit",true);
%%
%[text] ## Plot gene expression profile of cells ordered according to their pseudotime T.
r=corr(t,X','type','spearman'); % Calculate linear correlation between gene expression profile and T
[~,idxp]= maxk(r,4);  % Select top 4 positively correlated genes
[~,idxn]= mink(r,3);  % Select top 3 negatively correlated genes
selectedg=genelist([idxp idxn]);

% Plot expression profile of the 5 selected genes
try
    figure;
    gui.i_plot_pseudotimeseries(log1p(X),genelist,t,selectedg)
catch
end
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
% gui.i_plot_pseudotimeseries(log1p(X),genelist,t,selectedg)
%%
%[text] ## Trajectory analysis using TSCAN
%[text] Calculte pseudotime T
try
figure;
t=sc_trajectory(X,"type","tscan","plotit",true);

r=corr(t,X','type','spearman'); % Calculate linear correlation between gene expression profile and T
[~,idxp]= maxk(r,4);  % Select top 4 positively correlated genes
[~,idxn]= mink(r,3);  % Select top 3 negatively correlated genes
selectedg=genelist([idxp idxn]);
catch ME
    disp(ME.message);
end
% Plot expression profile of the 5 selected genes
try
figure;
gui.i_plot_pseudotimeseries(log1p(X),genelist,t,selectedg)
catch ME
    disp(ME.message);
end
%%
%[text] ## Construct single-cell gene regulatory network (scGRN)
%[text] ## Using principal component regression (PCNet) method
X50=X(1:50,:);
genelist50=genelist(1:50);
A=sc_grn(X50, 'pcrnet');

% Plot constructed network
%
A=A.*(abs(A)>quantile(abs(A(:)),0.9));
G=digraph(A,genelist50);
LWidths=abs(5*G.Edges.Weight/max(G.Edges.Weight));
LWidths(LWidths==0)=1e-5;
figure;
p=plot(G,'LineWidth',LWidths);
p.MarkerSize = 7;
p.Marker = 's';
p.NodeColor = 'r';
%%
%[text] ## Using GENIE3 method
X20=X(1:20,:);
genelist20=genelist(1:20);
A=run.ml_GENIE3(X20);

% Plot constructed network
%
A=A.*(abs(A)>quantile(abs(A(:)),0.9));
G=digraph(A,genelist20);
LWidths=abs(5*G.Edges.Weight/max(G.Edges.Weight));
LWidths(LWidths==0)=1e-5;
figure;
p=plot(G,'LineWidth',LWidths);
p.MarkerSize = 7;
p.Marker = 's';
p.NodeColor = 'r';
%%
%[text] ## The End

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
