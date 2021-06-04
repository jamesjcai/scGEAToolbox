function [T,Xsorted,genelistsorted]=sc_hvg(X,genelist,sortit,plotit,normit,ignorehigh)
% HVGs selection - This method uses the CV^2 on normalized count data to 
% analyze biological variation. 
%
% REF: https://www.nature.com/articles/nmeth.2645
% Input X: Brennecke et al. (2013) uses DESeq's method for normalization.
% 
% USAGE:
% >> [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% >> [T]=sc_hvg(X,genelist);


if nargin<2 || isempty(genelist)
    genelist=strcat("G",string(1:size(X,1)))';
end
if nargin<3, sortit=true; end
if nargin<4, plotit=false; end
if nargin<5, normit=true; end
if nargin<6, ignorehigh=true; end

if nargout>1, Xori=X; end

if normit
    [X]=norm_deseq(X);
    %[X]=norm_libsize(X);	
end
if any(isnan(X(:)))
    u=nanmean(X,2);
    vx=nanvar(X,0,2);
    cv2=vx./u.^2;
else
    u=mean(X,2);
    vx=var(X,0,2);
    % vx=sum(abs(X-mean(X,2)).^2,2)./(size(X,2)-1);
    cv2=vx./u.^2;
end

if issparse(u), u=full(u); end
if issparse(vx), vx=full(vx); end
if issparse(cv2), cv2=full(cv2); end

xi=1./u;
yi=cv2;

if ignorehigh
    yi=yi(xi>0.1);
    xi=xi(xi>0.1);
end
m=size(X,2);
df=m-1;

% b=glmfit(xi,yi,'gamma','link','identity');
% cv2fit=glmval(b,1./u,'identity');    % OR cv2fit=b(2)./u+b(1);
% if issparse(xi), xi=full(xi); end
% if issparse(yi), yi=full(yi); end
mdl=fitglm(xi,yi,'linear','Distribution','gamma','link','identity');
cv2fit=mdl.predict(1./u);
b=mdl.Coefficients.Estimate;


methodid=1;
switch methodid
    case 1
        % this code follows: https://github.com/tallulandrews/M3Drop/blob/master/R/Brennecke_implementation.R
        % and pages 7-9 https://media.nature.com/original/nature-assets/nmeth/journal/v10/n11/extref/nmeth.2645-S2.pdf
        minBiolDisp = 0.5.^2;
        cv2th = b(1) + minBiolDisp + b(1) * minBiolDisp;
        testDenom =(u*b(2) + u.^2*cv2th)/(1+cv2th/m);
        fitratio=vx./testDenom;
    case 2
        % this code follows https://github.com/MarioniLab/MNN2017/blob/a202f960f165816f22dec3b62ce1c7549b3ba8c1/Pancreas/findHighlyVariableGenes.R       
        % cv2fit=b(2)./u+b(1);
        fitratio=cv2./cv2fit;
end

pval=chi2cdf(fitratio*df,df,'upper');
% OR 1-chi2cdf(fitratio*df,df);
residualcv2=log(fitratio);   % log(cv2)-log(cv2fit);

% fdr=mafdr(pval,'BHFDR',true);
[~,~,~,fdr]=fdr_bh(pval);

T=table(genelist,u,cv2,residualcv2,fitratio,pval,fdr);
T.Properties.VariableNames(1)={'genes'};
i=~isnan(cv2);
T=T(i,:);
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
    figure;
    scatter(log(u),log(cv2));
    hold on
    % scatter(log(u(top100idx)),log(cv2(top100idx)),'x');
    plot(log(u),log(cv2fit),'.','markersize',10);    
    
    %[~,i]=sort(fitratio,'descend');
    %xi=u(i); yi=cv2(i); yifit=cv2fit(i);    
    %
%    scatter(log(xi),log(yi))
%    hold on
%    scatter(log(xi(1:100)),log(yi(1:100)),'x');
%    plot(log(xi),log(yifit),'.','markersize',10);   
%    plot(log(xi),log(yifit*chi2inv(0.975,df)./df),'.k');
%    plot(log(xi),log(yifit*chi2inv(0.025,df)./df),'.k');

    xlabel('Mean expression, log')
    ylabel('CV^2, log')
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
end