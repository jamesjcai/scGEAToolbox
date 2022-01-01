function scatter_meanlg_vs_varlg(X,genelist,dofit)
if nargin<2, genelist=[]; end
if nargin<3, dofit=false; end

Xn=log(X);
Xn(isnan(Xn)|isinf(Xn))=nan;
u=mean(Xn,2,'omitnan');
v=var(Xn,[],2,'omitnan');
scatter(u,v);

if dofit
    [xData, yData] = prepareCurveData( u, v );
    ft = fittype( 'a*x/(x^n+b)', 'independent', 'x', 'dependent', 'y' );
    fo = fitoptions( 'Method', 'NonlinearLeastSquares' );
    warning('off','curvefit:fit:noStartPoint');
    [fr] = fit( xData, yData, ft, fo );
    warning('on','curvefit:fit:noStartPoint');

    rangev=xlim();
    rangev(2)=max(xData);
    hold on
    i=rangev(1):0.05:rangev(2);
    plot(i,fr(i),'-r');
    %plot( fr, xData, yData );
end
grid on
xlabel('Mean(log(Expr))');
ylabel('Variance(log(Expr))');
if dofit
    legend({'Genes',sprintf('Nonlinear least squares\na*x/(x^n+b)')})
end
if ~isempty(genelist)
    dt = datacursormode;
    dt.UpdateFcn = {@i_myupdatefcn1,genelist};
end
