function [T, Xsorted_completed, gsorted_completed, ...
    xyz1] = sc_splinefit(X, genelist, sortit, plotit, removenan)
%SC_SPLINEFIT identify genes with a profile deviated from normal
%
% USAGE:
% >> [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% >> [X]=sc_norm(X,'type','libsize');
% >> [T]=sc_splinefit(X,genelist,true,true);

% if nargin<2 || isempty(genelist)
%     genelist=string(1:size(X,1))';
%     genelist=strcat("gene_",genelist);
% end

if nargin < 5, removenan = false; end
if nargin < 4, plotit = false; end
if nargin < 3, sortit = true; end
if nargin < 2 || isempty(genelist)
    genelist = string(1:size(X, 1)); 
end

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

% lgu=zscore(lgu);
% dropr=zscore(dropr);
% lgcv=zscore(lgcv);

%[~,i]=max(lgcv);
xyz = [lgu, lgcv, dropr];

% [~,j]=sort(pdist2(xyz,xyz(i,:)));
% xyz=xyz(j,:)';
% lgu=lgu(j);
% dropr=dropr(j);
% lgcv=lgcv(j);

%xyz=[lgu dropr lgcv]';

%s = cumsum([0;sqrt(diff(lgu(:)).^2 + diff(dropr(:)).^2 ...
%    + diff(lgcv(:)).^2)]);
s = cumsum([0; sqrt(diff(lgu(:)).^2+diff(lgcv(:)).^2 ...
    +diff(dropr(:)).^2)]);

pp1 = splinefit(s, xyz.', 15, 0.75);
xyz1 = ppval(pp1, s)';


[~, d] = dsearchn(xyz1, xyz);
fitmeanv=xyz1(:,1);
x=xyz(:,1); y=xyz(:,2);
d(x>max(fitmeanv))=d(x>max(fitmeanv))./100;
d(x<min(fitmeanv))=d(x<min(fitmeanv))./10;
d((y-xyz1(:, 2))<0)=d((y-xyz1(:, 2))<0)./100;


%D = pdist2(xyz, xyz1);
%d = min(D, [], 2);
dx = d(d <= quantile(d, 0.9));

distFit = fitdist([-dx; dx], 'Normal');
pval = normcdf(d, 0, distFit.sigma, 'upper');
[~, ~, ~, fdr] = pkg.fdr_bh(pval);

if ~isempty(gsorted)
    genes = gsorted;
    T = table(genes, lgu, lgcv, dropr, d, pval, fdr);
else
    T = table(lgu, lgcv, dropr, d, pval, fdr);
end
% 'variablenames',{'Genes','Log10_Mean','Dropout_Rate','Log10_CV','Deviation_3DFeature'});

T.d(T.dropr > (1 - 0.05)) = 0; % ignore genes with dropout rate > 0.95
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
