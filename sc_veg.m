function [T,mdl]=sc_veg(X,genelist,sortit,plotit,donorm)
%SC_VEG 
%
% REF: https://github.com/hillas/scVEGs
% 
% USAGE:
% >> [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% >> [X]=sc_norm(X,'type','libsize');
% >> [T]=sc_veg(X,genelist,true,true,false);

if nargin<5, donorm=false; end
if nargin<4, plotit=false; end
if nargin<3, sortit=true; end

if donorm
    [X]=sc_norm(X,'type','libsize');
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

fitm = locfit(lgxdata,ydata);

% figure;
% lfplot(fitm);
% hold on
% ysmooth = malowess(lgxdata,ydata);
% [~,idx]=sort(lgxdata);
% plot(sort(lgxdata),ysmooth(idx),'-r');

xSeq=min(log10(xdata)):0.005:max(log10(xdata));
gapNum=zeros(length(xSeq),1);


for i=1:length(xSeq)
    gapNum(i)=sum((log10(xdata) >= xSeq(i) - 0.05) & (log10(xdata) < xSeq(i) + 0.05));
end
cdx=find(gapNum > m*0.005);

  xSeq = 10.^ xSeq;
%  ySeq = predict(fitm,log10(xSeq));
  ySeq=zeros(1,length(xSeq));
  for k=1:length(ySeq)
      ySeq(k)=predict(fitm,log10(xSeq(k)));
  end
  
yDiff=diff(ySeq);
ix=find(yDiff>0 & log10(xSeq(2:end))>0);
if isempty(ix), ix=length(ySeq)-1; end

xSeq_all=10.^(min(log10(xdata)):0.001:max(log10(xdata)));
xSeq=xSeq((cdx(1):ix(1))+1);
ySeq=ySeq((cdx(1):ix(1))+1);


% Fit model to data.
[xData, yData] = prepareCurveData( xSeq, ySeq );
ft = fittype( '0.5*log10(b/x+a)', 'independent', 'x', 'dependent', 'y' );
fo = fitoptions( 'Method', 'NonlinearLeastSquares' );
fo.StartPoint = [0 1];
% opts.Display = 'Off';
[mdl] = fit( xData, yData, ft, fo );
ydataFit = mdl(xSeq_all)';


% Calculate CV difference
  logX = log10(xdata);
  logXseq = log10(xSeq_all);
  cvDist =zeros(length(xdata),1);
  
  for i=1:length(logX)
    cx = find(logXseq >= logX(i) - 0.2 & logXseq < logX(i) + 0.2);
    tmp = sqrt((logXseq(cx) - logX(i)).^2 + (ydataFit(cx) - ydata(i)).^2);
    [~,tx]=min(tmp);
    tx=tx(1);
    if(logXseq(cx(tx)) > logX(i)) 
      if(ydataFit(cx(tx)) > ydata(i))
        cvDist(i) = -tmp(tx);
      else
        cvDist(i) = tmp(tx);
      end
      cvDist(i) = -tmp(tx);
    elseif (logXseq(cx(tx)) <= logX(i))
      if(ydataFit(cx(tx)) < ydata(i))
        cvDist(i) = tmp(tx);
      else
        cvDist(i) =-tmp(tx);
      end
    end
  end
cvDist = log2(10.^cvDist);



[f,xi] = ksdensity(cvDist,'NumPoints',512);

[~,idx]=max(f);
distMid=xi(idx);
dist2 = cvDist - distMid;
tmpDist = [dist2(dist2 <= 0); abs(dist2(dist2 < 0))] + distMid;

distFit = fitdist(tmpDist,'Normal');
pval=normcdf(cvDist,distFit.mu,distFit.sigma,'upper');

[~,~,~,fdr]=fdr_bh(pval);

residualcv=cvDist;
T=table(genelist,u,cv,residualcv,pval,fdr);
T.Properties.VariableNames(1)={'genes'};
if sortit
    T=sortrows(T,'residualcv','descend');
end

if plotit
    % [~,top100idx]=maxk(fitratio,100);
    scatter(log10(u),log10(cv));
    hold on
    % scatter(log(u(top100idx)),log(cv2(top100idx)),'x');
    % plot(log(u),log(cv2fit),'.','markersize',10);    

  %plot(log10(xSeq), ySeq,'.','markersize',10);
  plot(log10(xSeq_all), ydataFit,'.','markersize',10);
    
    xlabel('Mean expression, log')
    ylabel('CV, log')
    if ~isempty(genelist)
        dt=datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn3,genelist,X};
    end
    hold off
end

end