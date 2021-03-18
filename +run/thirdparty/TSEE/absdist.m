function absd = absdist(X_time)
N = size(X_time,1);
absd = abs(bsxfun(@minus, X_time,X_time'));
