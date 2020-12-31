function [T]=sc_lcod(X,genelist,sortit,plotit,donorm)
% Lowess coefficient of dispersion (LCOD) analysis
% To select the variable genes we first fitted a Lowess curve to the log2 of the mean
% versus the log2 of the standard deviation of expression among cells and then
% calculated for each gene the distance from this curve. The distribution of these
% distances is approximately normal distributed; for this reason we used a z-score to
% select the most and the least variable genes (>1.5 and <-1.5 respectively)
%
% REF: https://www.cell.com/cms/10.1016/j.celrep.2015.12.089/attachment/15ca54cd-72cc-44c4-a421-437ed5f4d733/mmc1
% Input X: 
% 
% USAGE:
% >> [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% >> [X]=sc_norm(X,'type','libsize');
% >> [T]=sc_lcod(X,genelist,true,true,false);

if nargin<5, donorm=false; end
if nargin<4, plotit=false; end
if nargin<3, sortit=true; end

if donorm
    [X]=norm_libsize(X);
end
u=nanmean(X,2);
sd=nanstd(X,0,2);


xi=log2(u);
yi=log2(sd);

% m=size(X,2);
% df=m-1;

% [dataout]=lowess([xi yi],0.25);
[yi_smoothed]=fLOESS([xi yi],15);
res_yi=yi-yi_smoothed;
% i=abs(zscore(res_yi))>1.5;
distFit = fitdist(res_yi,'Normal');
pval=normcdf(yi,distFit.mu,distFit.sigma,'upper');

[~,~,~,fdr]=fdr_bh(pval);

%%
residual_sd=res_yi;
T=table(genelist,u,sd,residual_sd,pval,fdr);
T.Properties.VariableNames(1)={'genes'};
if sortit
    T=sortrows(T,'residual_sd','descend');
end

if plotit
    % [~,top100idx]=maxk(fitratio,100);
    scatter(xi,yi);
    hold on
    % scatter(log(u(top100idx)),log(cv2(top100idx)),'x');
    % plot(log(u),log(cv2fit),'.','markersize',10);    

  plot(xi, yi_smoothed,'.','markersize',10);
      
    xlabel('Mean, log2')
    ylabel('SD, log2')
    if ~isempty(genelist)
        dt=datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn3,genelist,X};
    end
    hold off
end