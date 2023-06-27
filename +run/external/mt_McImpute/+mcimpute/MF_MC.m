function [X, U, V] = MF_MC(y,M,sizeX,rankr)

% Matrix Completion via Factorization
% min nuclear-norm(X) subject to ||y - M(X)||_2<err

% Inputs
% X - matrix to be estimated
% M - masking operator, applied to vectorized form of X
% y - sampled entries



err = 1e-6;
x = rand(prod(sizeX),1);
outsweep = 10;
insweep = 20;
tol = 1e-4;    
decfac = 0.9;
alpha = 1.1;

[U1,S1,V1] = svd(reshape(x,sizeX));
V = S1(1:rankr,1:rankr)*V1(:,1:rankr)';

for out = 1:outsweep
        x = x + (1/alpha)*M(y - M(x,1),2);
        X = reshape(x,sizeX);
    for ins = 1:insweep
        U = mrdivide(X, V);
        for j = 1:rankr
            U(:,j) = U(:,j)/norm(U(:,j));
        end
        V = mldivide(U,X);
    end
        X = U*V;
        X(X<0)=0;
        x = X(:);
end
