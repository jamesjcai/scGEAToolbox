function [X]  = IST_MC(y,M,sizeX,rankX)

% Matrix Completion via Iterated Soft Thresholding
% min nuclear-norm(X) subject to ||y - M(X)||_2<err

% Inputs
% X - matrix to be estimated
% M - masking operator, applied to vectorized form of X
% y - sampled entries
% err - norm of the mismatch (default 0)
% x_initial - intial estimate of the vectorized form of X (defalut 0)
% normfac - eigenvalue of (M'M) (default should be 1 for masking operator)
% insweep - maximum number of internal sweeps for solving ||y - M(X)||_2 + lambda nuclear-norm(X) (default 200)
% tol - tolerance (default 1e-4)
% decfac - decrease factor for cooling lambda

    err = 1e-12;
    x_initial = zeros(prod(sizeX),1); % leave it as it is
    normfac = 1; % leave it as it is
    insweep = 20;
    tol = 1e-4;    
    decfac = 0.7;%earlier it was 0.9

alpha =1.1*normfac;
x = x_initial;
lambdaInit = decfac*max(abs(M(y,2))); lambda = lambdaInit;
f_current = norm(y-M(x,1)) + lambda*norm(x,1);

while lambda >lambdaInit*tol
    lambda
    for ins = 1:insweep        
        f_previous = f_current;
        b = x + (1/alpha)*M(y - M(x,1),2); % Landweber Iteration
        B = reshape(b,sizeX);
        [U,S,V] =svd(B,'econ'); %lansvd(B,rankX,'L'); %  uses more efficient lansvd from PROPACK
        s = SoftTh(diag(S),lambda/(2*alpha));
        S = diag(s);
        X = U*S*V';
        X(X<0)=0;
        x = X(:);
        f_current = norm(y-M(x,1)) + lambda*norm(x,1);
        if norm(f_current-f_previous)/norm(f_current + f_previous)<tol
            break;
        end
    end
    if norm(y-M(x,1))<err
        break;
    end
    
    lambda = decfac*lambda;
end
end

function  z = SoftTh(s,thld)
    z = sign(s).*max(0,abs(s)-thld); 
end

