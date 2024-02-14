function [T, Xsorted_completed, gsorted_completed, ...
    xyz1] = sc_splinefit(X, genelist, sortit, plotit, removenan)
%{
    SC_SPLINEFIT identify genes with a profile deviated from normal
    INPUTS
        X ----------> Count matrix
        genelist ---> gene list for count matrix X 
        sortit -----> boolean to sort or not the gene statistics (it's required)
        plotit -----> boolean to plot cubic spline fit
        removenan --> removing noisy and undesired NaN data
    OUTPUT
        T ----------> table containing log mean (lgu), log coef_variance (lgcv),
                      drop rate (dropr), distance spline-gene (d), pvalue
                      (pval) and false discovery rate (fdr)
        Xsorted ----> Sorted Xcount matrix according to T order (completed)
        gsorted ----> Sorted genelist according to T order (completed)
        xyz1 -------> Spline fit (cubic spline) of noisy data evaluated
                      into s-grid-size (lgu size)
    USAGE
        [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
        [X]=sc_norm(X,'type','libsize');
        [T]=sc_splinefit(X,genelist,true,true);
%}

if nargin < 5, removenan = true; end
if nargin < 4, plotit = false; end
if nargin < 3, sortit = true; end
if nargin < 2, genelist = string(1:size(X, 1)); end

% Computing log-mean, log-coef_variance and drop out rate 
[lgu, dropr, lgcv, gsorted, Xsorted, ...
    removedgidx, removedT] = sc_genestat(X, genelist, sortit, removenan);

if removenan && ~isempty(removedgidx)
    gsorted_completed = [gsorted; genelist(removedgidx)];
    Xsorted_completed = [Xsorted; X(removedgidx,:)];
    assert(isequal(size(Xsorted_completed),size(X)),'SC_SPLINEFIT')
    assert(length(gsorted_completed) == length(genelist),'SC_SPLINEFIT')
else
    gsorted_completed = gsorted;
    Xsorted_completed = Xsorted;
end

% xyz contains log-mean, log-coef_variance and drop out rate 
xyz = [lgu, lgcv, dropr];

% s is what? statistics curve?
% Cumulative sum( sqrt( D1^2 + D2^2 + D3^2) )
%  Di = [Xi(2)-Xi(1) Xi(3)-Xi(2) ... Xi(m)-Xi(m-1)]
%  Xi can be lgu_i, lgcv_i, dropr_i
s = cumsum( [ 0; sqrt( diff( lgu(:)  ).^2 + ...
                       diff( lgcv(:) ).^2 + ...
                       diff( dropr(:)).^2) ] );

% Construct polinomial points of cubic spline of noisy data
pp1 = splinefit(s, xyz.', 15, 0.75);
% Polynomial points evaluated from splinefit into s-grid-size
xyz1 = ppval(pp1, s)';

% search for nearest points
[~, d] = dsearchn(xyz1, xyz);

% passing the polynomial fit from gene statistics
fitmeanv = xyz1(:,1);

% x and y values from lgu and lgcv 
x = xyz(:,1); y = xyz(:,2);

% if x is far from mean curve, then scale?
d(x>max(fitmeanv)) = d(x>max(fitmeanv))./100;
d(x<min(fitmeanv)) = d(x<min(fitmeanv))./10;
d((y-xyz1(:, 2))<0) = d((y-xyz1(:, 2))<0)./100;

% Obtain distances that are up to a cumulative probability up to 90%
dx = d(d <= quantile(d, 0.9));

% Normal distribution fit for distances in dx
distFit = fitdist([-dx; dx], 'Normal');

% Cumulative distribution function (cdf) according to
% a standard normal distribution in d from dx sigma
pval = normcdf(d, 0, distFit.sigma, 'upper');
[~, ~, ~, fdr] = pkg.fdr_bh(pval);

if ~isempty(gsorted)
    genes = gsorted;
    T = table(genes, lgu, lgcv, dropr, d, pval, fdr);
else
    T = table(lgu, lgcv, dropr, d, pval, fdr);
end

% 'variablenames',{'Genes','Log10_Mean','Dropout_Rate','Log10_CV','Deviation_3DFeature'});

% ignore genes with dropout rate > 0.95
T.d(T.dropr > (1 - 0.05)) = 0; 
% disp('NOTE: Genes with dropout rate > 0.95 are excluded.');

if ~isempty(removedT) && istable(removedT)
    removedT.Properties.VariableNames = T.Properties.VariableNames;
    T = [T; removedT];
end

if sortit
    [T,idx] = sortrows(T, 'd', 'descend');
    gsorted_completed = gsorted_completed(idx);
    Xsorted_completed = Xsorted_completed(idx,:);
end

if length(gsorted_completed) ~= length(genelist)
    error('Output GENES are less than input GENES (some GENES are removed).');
end

if plotit
    figure;
    scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 'filled', 'MarkerFaceAlpha', .1);
    hold on
    plot3(xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), '-', 'linewidth', 4);
    xlabel('Mean, log');
    ylabel('CV, log');
    zlabel('Dropout rate (% of zeros)');

    if ~isempty(gsorted)
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn3, gsorted, Xsorted};
    end
    hold off
end

end
