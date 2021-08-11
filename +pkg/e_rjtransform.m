function [Y]=e_rjtransform(X)
% RJclust method transformation
    [m,n]=size(X);
    if m~=n, error('X should be a square matrix.'); end
    Z1=tril(X,-1);
    Z2=triu(X,1);
    Y=[Z1(:,1:end-1)+Z2(:,2:end) diag(X)];
end