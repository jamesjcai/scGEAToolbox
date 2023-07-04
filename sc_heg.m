function [T,Xsorted,genelistsorted]=sc_heg(X,genelist,~,plotit,normit,~)
% HEGs - highly expressed genes

if nargin<2 || isempty(genelist)
    genelist=strcat("G",string(1:size(X,1)))';
end
if nargin<3, sortit=true; end
if nargin<4, plotit=true; end
if nargin<5, normit=true; end
if nargin<6, ignorehigh=true; end

if nargout>1, Xori=X; end

if normit
    [X]=pkg.norm_deseq(X);
end

u=mean(X,2,'omitnan');
cv=std(X,0,2,'omitnan')./u;

m=size(X,1);
ydata=log10(cv);            
xdata=u(~isnan(ydata));

ydata=ydata(~isnan(ydata));
lgxdata=log10(xdata);


pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'+run','thirdparty','locfit','m');
if ~(ismcc || isdeployed), addpath(pth); end
pth=fullfile(pw1,'+run','thirdparty','locfit','mex');
if ~(ismcc || isdeployed), addpath(pth); end
pth=fullfile(pw1,'+run','thirdparty','locfit','source');
if ~(ismcc || isdeployed), addpath(pth); end

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
if plotit
figure;
    scatter(ydata,lgxdata)
    hold on
    scatter(ydata,fit_lgxdata)
end

end

