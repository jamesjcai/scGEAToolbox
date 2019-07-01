function scatter_mean_vs_cv(X,genelist,dofit)
if nargin<2, genelist=[]; end
if nargin<3, dofit=false; end

u=nanmean(X,2);
cv=nanstd(X,[],2)./u;

loglog(u,cv,'o');
grid on
xlabel('u, log10');
ylabel('CV, log10');
% labels following https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4441768/
% xlabel('avg. counts (UMI filtered), log10');
% ylabel('inter-cell CV, log10');

if dofit
    lgcv=log10(cv);
    [xData, yData] = prepareCurveData( u, lgcv );
    ft = fittype( '0.5*log10(b/x+a)', 'independent', 'x', 'dependent', 'y' );
    fo = fitoptions( 'Method', 'NonlinearLeastSquares' );
    warning('off','curvefit:fit:noStartPoint');
    [fr] = fit( xData, yData, ft, fo );
    warning('on','curvefit:fit:noStartPoint');
    % rangev=log10([min(xData) max(xData)]);
    rangev=log10(xlim());
    rangev(2)=log10(quantile(xData,0.99));
    ab=coeffvalues(fr);
    hold on

    i=10^rangev(1):0.04:10^rangev(2);

    j=10.^-(0.5*log10(i));
    plot(i,j,'-rv');
    j=10.^(0.5*log10(ab(2)./i+ab(1)));
    % plot(i,j,'-gs');
    % plot(log10(sort(xData)),fr(sort(xData)),'g-');

    %i=logspace(rangev(1),rangev(2));
    j=10.^(0.5*log10(ab(2)./i+ab(1)));
    plot(i,j,'-gs');
    legend({'Genes','Poisson distribution',...
        sprintf('Nonlinear least squares\n(a=%.3f, b=%.3f)',...
        ab(1),ab(2))},'location','southwest')
end

if ~isempty(genelist)
    dt = datacursormode;
    dt.UpdateFcn = {@i_myupdatefcn1,genelist};
end

