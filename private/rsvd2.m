function [U,S,V]=rsvd2(X,r,q,p)
    % https://youtu.be/AFyfEqh6g3I?t=353
    if nargin<4, p=5; end
    if nargin<3, q=5; end
    if nargin<2, r=50; end
    ny=size(X,2);
    P=randn(ny,r+p);
    Z=X*P;
    for k=1:q
        Z=X*(X'*Z);   % power iterations
    end
    [Q,~]=qr(Z,0);

    Y=Q'*X;
    [UY,S,V]=svd(Y,'econ');
    U=Q*UY;
end
