function [T,Xsorted,genelistsorted]=sc_heg(X,genelist,sortit,plotit,normit,ignorehigh)
% HEGs - highly expressed genes
%
%
% REF: https://www.nature.com/articles/nmeth.2645
% Input X: Brennecke et al. (2013) uses DESeq's method for normalization.
% 
% USAGE:
% >> [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% >> [X]=sc_norm(X,'type','deseq');
% >> [T]=sc_hvg(X,genelist,true,true);


if nargin<2 || isempty(genelist)
    genelist=strcat("G",string(1:size(X,1)))';
end
if nargin<3, sortit=true; end
if nargin<4, plotit=true; end
if nargin<5, normit=true; end
if nargin<6, ignorehigh=true; end

if nargout>1, Xori=X; end

if normit
    [X]=norm_deseq(X);
end

u=nanmean(X,2);
cv=nanstd(X,0,2)./u;

m=size(X,1);
ydata=log10(cv);            
xdata=u(~isnan(ydata));

ydata=ydata(~isnan(ydata));
lgxdata=log10(xdata);


pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty/locfit/m');
addpath(pth);
pth=fullfile(pw1,'thirdparty/locfit/mex');
addpath(pth);
pth=fullfile(pw1,'thirdparty/locfit/source');
addpath(pth);

%xi=ydata;
%yi=lgxdata;
fit_lgxdata=zeros(size(lgxdata));
fitm = locfit(ydata,lgxdata);    % log10(cv) vs log10(u)
  for k=1:length(ydata)
      fit_lgxdata(k)=predict(fitm,ydata(k));
  end
Xsorted=[];
genelistsorted=[];
residualu=fit_lgxdata-lgxdata;   % fitted_u-u   
T=table(genelist,u,cv,residualu);

T=T(T.u>quantile(T.u,0.4),:);

figure;
scatter(ydata,lgxdata)
hold on
scatter(ydata,fit_lgxdata)
return;
%{
u=nanmean(X,2);
cv2=nanvar(X,0,2)./u.^2;

% xi=log(cv2);
% yi=log(u);
% mdl = fitlm(xi,yi);
% scatter(xi,mdl.predict)




%xi=1./u;
%yi=cv2;
xi=cv2;
yi=u;

if ignorehigh
   % yi=yi(xi>0.02);
   % xi=xi(xi>0.02);
   yi=yi(xi<50);
   xi=xi(xi<50);
   
end
m=size(X,2);
df=m-1;

% b=glmfit(xi,yi,'gamma','link','identity');
% cv2fit=glmval(b,1./u,'identity');    % OR cv2fit=b(2)./u+b(1);

mdl=fitglm(xi,yi,'linear','Distribution','gamma','link','identity');
%cv2fit=mdl.predict(1./u);
u_fit=mdl.predict(cv2);
fitratio=u./u_fit;
% b=mdl.Coefficients.Estimate;

pval=chi2cdf(fitratio*df,df,'upper');
% OR 1-chi2cdf(fitratio*df,df);
residualu=log(fitratio);   % log(cv2)-log(cv2fit);

% fdr=mafdr(pval,'BHFDR',true);
[~,~,~,fdr]=fdr_bh(pval);

T=table(genelist,u,cv2,residualu,fitratio,pval,fdr);
T.Properties.VariableNames(1)={'genes'};
if sortit
    [T,idx]=sortrows(T,'fitratio','descend');    
    if nargout>1
        Xsorted=Xori(idx,:);        
    end
    if nargout>2
        genelistsorted=T.genes;
    end
else
    if nargout>1
        error('SORTIT was not required.');
    end
end

if plotit
    % [~,top100idx]=maxk(fitratio,100);
    scatter(log(cv2),log(u));
    hold on
    % scatter(log(u(top100idx)),log(cv2(top100idx)),'x');
    plot(log(cv2),log(u_fit),'.','markersize',10);    
    
    %[~,i]=sort(fitratio,'descend');
    %xi=u(i); yi=cv2(i); yifit=cv2fit(i);    
    %
%    scatter(log(xi),log(yi))
%    hold on
%    scatter(log(xi(1:100)),log(yi(1:100)),'x');
%    plot(log(xi),log(yifit),'.','markersize',10);   
%    plot(log(xi),log(yifit*chi2inv(0.975,df)./df),'.k');
%    plot(log(xi),log(yifit*chi2inv(0.025,df)./df),'.k');

xlabel('CV^2, log')    
ylabel('Mean expression, log')
    
    if ~isempty(genelist)
        dt=datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn3,genelist,X};
    end
    hold off
end

% Highly variable genes (HVG) is based on the assumption that genes with 
% high variance relative to their mean expression are due to biological 
% effects rather than just technical noise. The method seeks to identify 
% genes that have a higher variability than expected by considering the 
% relationship between variance and mean expression. This relationship is 
% difficult to fit, and in practice genes are ranked by their distance 
% from a moving median (Kolodziejczyk et al., 2015) or another statistic 
% derived from variance is used, e.g. the squared coefficient of variation
% (Brennecke et al. (2013)).
%}