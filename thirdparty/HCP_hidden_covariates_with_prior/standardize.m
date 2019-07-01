function [X,m,st] = standardize(X)

[n p] = size(X);

%X = X - nanmean(X,2)*ones(1,p);
X=X-nanmean(X,2);

m = mean(X,2);
X(isnan(X)) = 0;
st = sqrt(sum(X.^2,2));
X = X./sqrt(sum(X.^2,2)*ones(1,p));
