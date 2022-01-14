function [t]=TSCAN(X,varargin)
%TSCAN - pseudotime analysis

% ref: https://academic.oup.com/nar/article/44/13/e117/2457590
% load example_data\example10xdata.mat
% [X,genelist]=sc_selectg(X,genelist,7,8);
% [X]=sc_selectc(X);
% [t_pseudotime]=sc_tscan(X);

p = inputParser;
addRequired(p,'X',@isnumeric);
addOptional(p,'do_geneculst',true,@islogical);
addOptional(p,'do_reduce',true,@islogical);
addOptional(p,'plotit',true,@islogical);
parse(p,X,varargin{:})

do_geneculst=p.Results.do_geneculst;
do_reduce=p.Results.do_reduce;
plotit=p.Results.plotit;

if do_geneculst
    % in order to alleviate the effect of drop-out events (20) on 
    % the subsequent analyses, genes with similar expression patterns are 
    % grouped into clusters by hierarchical clustering (using Euclidean 
    % distance and complete linkage). The number of clusters is set to 
    % be 5% of the total number of genes with non-zero expression. For each 
    % cluster and each cell, the expression measurements of all genes in 
    % the cluster are averaged to produce a cluster-level expression which 
    % will be used for subsequent MST construction. 
    n=round(0.05*size(X,1));
    Y=pdist(X);
    Z=linkage(Y);
    clu=cluster(Z,'maxclust',n);
    X=grpstats(X,clu,@(x) mean(x,1));
    % Xn=[];
    % for k=1:27
    %    Xn=[Xn;mean(X(idx==k,:),1)];
    % end    
end



if do_reduce
    % Briefly, Ei from all cells are organized into a H × N matrix E?. Each row
    % corresponds to a gene cluster. The matrix is standardized such that 
    % expression values within each row have zero mean and unit standard deviation
    pcadim=5;
    Xn=zscore(X,0,2);
    [~,X] = pca(Xn','NumComponents',pcadim);
end


% z A matrix whose [i,k]th entry is the probability that observation i in the test data
% belongs to the kth class. https://github.com/zji90/TSCAN/blob/master/R/exprmclust.R
% G The optimal number of mixture components
% There are other options you can use to help select the appropriate number of components for a Gaussian mixture model. For example,
% Compare multiple models with varying numbers of components using information criteria, e.g., AIC or BIC.
% Estimate the number of clusters using evalclusters, which supports, the Calinski-Harabasz criterion and the gap statistic, or other criteria.
% res <- suppressWarnings(Mclust(pcareduceres, G = clusternum, modelNames = modelNames)
% https://www.mathworks.com/help/stats/clustering-using-gaussian-mixture-models.html

warning off
vaic=zeros(1,8);
% vbic=zeros(1,8);
rng(271)
for k=1:8
    try
        % gm=fitgmdist(X,k+1,'CovarianceType','full');
        gm=fitgmdist(X,k+1,'RegularizationValue',0.1);
        %vbic(k)=gm.BIC;
        vaic(k)=gm.AIC;
    catch
        %vbic(k)=nan;
        vaic(k)=nan;
    end
end
[~,idx1]=min(vaic);
% [~,idx2]=min(vbic);
% clunum=min([idx1 idx1])+1;
clunum=idx1+1;
% clunum=fun_num_cluster(X');
warning on


gmfit=fitgmdist(X,clunum,'CovarianceType','full');
clusterid = cluster(gmfit,X);
clucenter=grpstats(X,clusterid,@mean);

txtc=strings(size(clucenter,1),1);
for k=1:size(clucenter,1)
    txtc(k)=string(sprintf('Clu%d',k));
end

G=graph(squareform(pdist(clucenter)),txtc);
T=minspantree(G);
D=distances(T);
[i,j]=find(D==max(D(:)));
clupath=shortestpath(T,i(1),j(1));


%%
clupath=[clupath clupath(1)];
tt=nan(size(X,1),2);
for k=1:length(clupath)-1
    i=clupath(k);
    j=clupath(k+1);
    idx=clusterid==i;
    
    c1=clucenter(i,1:3);
    c2=clucenter(j,1:3);
    
    x1=X(idx,1:3);
    tt(idx,1)=k;
    
    difvec=c2-c1;
    difv=difvec/norm(difvec);
    [~,idxv]=sort(difv*x1');
    tt(idx,2)=idxv;
    
%     x1=x1(idx,:);
%     for k=1:size(x1,1)
%         text(x1(k,1),x1(k,2),sprintf('%d',k));
%     end
end
[~,t]=sortrows(tt,[1 2]);

if plotit    
    subplot(2,2,2)
    p = plot(G); % ,'EdgeLabel',G.Edges.Weight);    
    highlight(p,T,'EdgeColor','r','LineWidth',4.5)
    title('Graph of cell cluster centers')
    
    subplot(2,2,3)
    plot(T);
    title('MST')

    subplot(2,2,4)
    p=plot(T);
    highlight(p,clupath,'EdgeColor','r','LineWidth',4.5)
    title('Main path in MST')
    
    subplot(2,2,1)
    gui.i_gscatter3(X(:,1:3),clusterid);
    hold on
    for k=1:clunum
        plot3(clucenter(k,1),clucenter(k,2),clucenter(k,3),'+','markersize',20)
        text(clucenter(k,1),clucenter(k,2),clucenter(k,3),sprintf('%d',k),'fontsize',30)
    end
    title('Cell clusters')
end

end

%{
c1=clucenter(1,1:2);
c2=clucenter(2,1:2);
x1=X(clusterid==1,1:2);
difvec=c2-c1;
difv=difvec/norm(difvec);
[~,idx]=sort(difv*x1');
x1=x1(idx,:);
for k=1:size(x1,1)
    text(x1(k,1),x1(k,2),sprintf('%d',k));
end
%}
