%function [p]=deseq2a(X,Y,g)

tLocal1 = nbintest(X, Y, 'VarianceLink', 'LocalRegression');
p1 = tLocal1.pValue;
sfx1 = tLocal1.SizeFactors{1};
sfy1 = tLocal1.SizeFactors{2};

geom1 = geomean(matrix1);
geom1(geom1 == 0) = 1;
norm1 = median(bsxfun(@rdivide, matrix1, geom1), 2);

geom1 = geomean([matrix1; matrix2]);
%geom1(geom1==0)=1;
norm1 = median(bsxfun(@rdivide, [matrix1; matrix2], geom1), 2);
%matrixN1 = bsxfun(@rdivide,matrix1',norm1')';


[Xn, sfx] = norm_deseq(X);
[Yn, sfy] = norm_deseq(Y);
tLocal = nbintest(X, Y, 'VarianceLink', 'LocalRegression', ...
    'SizeFactor', {sfx', sfy'});
p = tLocal.pValue;
%end