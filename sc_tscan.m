% First, genes with zero read count in all samples are excluded. 
% load example_data\example10xdata.mat
% [X,genelist]=sc_selectg(X,genelist,8,7);
% [X]=sc_selectc(X);
load example_data\tscan_lpsdata.mat

% Second, in order to alleviate the effect of drop-out events (20) on 
% the subsequent analyses, genes with similar expression patterns are 
% grouped into clusters by hierarchical clustering (using Euclidean 
% distance and complete linkage). The number of clusters is set to 
% be 5% of the total number of genes with non-zero expression. For each 
% cluster and each cell, the expression measurements of all genes in 
% the cluster are averaged to produce a cluster-level expression which 
% will be used for subsequent MST construction. 
do_geneculst=true;
if do_geneculst
    n=round(0.05*size(X,1));
    Y=pdist(X);
    Z=linkage(Y);
    clu=cluster(Z,'maxclust',n);
    X=grpstats(X,clu,@(x) mean(x,1));
end
% Xn=[];
% for k=1:27
%    Xn=[Xn;mean(X(idx==k,:),1)];
% end

do_reduce=true;
if do_reduce
    % Briefly, Ei from all cells are organized into a H × N matrix E?. Each row
    % corresponds to a gene cluster. The matrix is standardized such that 
    % expression values within each row have zero mean and unit standard deviation
    [~,score,latent,~,explained] = pca(X');
    % figure()
    % pareto(explained)
    % xlabel('Principal Component')
    % ylabel('Variance Explained (%)')
    pcadim=5;
    Xn=zscore(X,0,2);
    [coeff,X] = pca(Xn','NumComponents',pcadim);
end

% z A matrix whose [i,k]th entry is the probability that observation i in the test data
% belongs to the kth class. https://github.com/zji90/TSCAN/blob/master/R/exprmclust.R
% G The optimal number of mixture components
%%

% There are other options you can use to help select the appropriate number of components for a Gaussian mixture model. For example,
%
% Compare multiple models with varying numbers of components using information criteria, e.g., AIC or BIC.
% Estimate the number of clusters using evalclusters, which supports, the Calinski-Harabasz criterion and the gap statistic, or other criteria.

vaic=zeros(1,8);
vbic=zeros(1,8);
for k=1:8
    try
        gm=fitgmdist(X,k+1,'CovarianceType','full');
        vbic(k)=gm.BIC;
        vaic(k)=gm.AIC;
    catch
        vbic(k)=nan;
        vaic(k)=nan;        
    end
end
[~,idx1]=min(vaic);
[~,idx2]=min(vbic);
clunum=min([idx1 idx2])+1;

clunum
clunum=3;

gmfit=fitgmdist(X,clunum,'CovarianceType','full');
% https://www.mathworks.com/help/stats/clustering-using-gaussian-mixture-models.html
clusterid = cluster(gmfit,X);
figure; i_myscatter(X(:,1:2),clusterid)
clucenter=grpstats(X,clusterid,@mean);
G=graph(squareform(pdist(clucenter)));
% https://academic.oup.com/nar/article/44/13/e117/2457590
T = minspantree(G);
%%
figure;
p = plot(G,'EdgeLabel',G.Edges.Weight);
highlight(p,T)
%%
figure;
i_myscatter(X(:,1:2),clusterid)
hold on
for k=1:clunum
    plot(clucenter(k,1),clucenter(k,2),'+','markersize',20)
    text(clucenter(k,1),clucenter(k,2),sprintf('%d',k),'fontsize',30)
end

