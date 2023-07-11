function [Z,B,U,o] = hidden_covariates_model(F,Y,k,lambda,lambda2,lambda3,iter)
%
%
% function [Z,B,U,o,error,error1,error2,dz,db,du] = hidden_covariate_linear(F,Y,k,lambda,lambda2,lambda3,iter);
% input:
%      F: a matrix nxd of known covariates, where n is the number of
%      subjects and d is the number of known covariates. *must be standardize (columns have 0 mean and constant SS).
%      Y: a matrix of nxg of expression data (must be standardized (columns
%      scaled to have constant SS and mean 0). ** use standardize function to standardize F and Y.
%      k: number of inferred hidden components (k is an integer)
%      lambda, lambda2, lambda3 are model parameters
%      (optional) iter: number of iterations (default = 100);
%
%      note: k>0, lambda>0, lambda2>0, lambda3>0 must be set by the user based on the data at hand. one can set these values
%      using cross-validation, by evaluating the "performance" of the  resulting residual data on a desired task.
%      typically, if lambda>5, then hidden factors match the known covariates closely. 
%
% objective:
%
% this function solves the following problem:
% argmin_{Z,B,U}   ||Y-Z*B||_2 + \lambda*||Z-F*U||_2 + \lambda2*||B||_2 + \lambda_3||U||_2
%
% output:
%      Z: matrix of hidden components, dimensionality: nxk
%      B: matrix of effects of hidden components, dimensionality: kxg
%      o: value of objective function on consecutive iterations.
%
% to use the residual data: Residual = Y - Z*B
%

if nargin < 7
   iter = 100;
end

% convergence criteria
tol = 1e-6;

U = zeros(size(F,2),k);
Z = zeros(size(F,1),k);
B = rand(size(Z,2),size(Y,2));

[n1,~] = size(F);
[n2,~] = size(Y);

if n1~=n2
    error('number of rows in F and Y must agree');
end

if (k<1 || lambda<1e-6 || lambda2<1e-6 || lambda3<1e-6 )
    error('lambda, lambda2, lambda3 must be positive and/or k must be an integer');
end



for ii = 1:iter
    
    o(ii) = sum(sum((Y-Z*B).^2)) + sum(sum((Z-F*U).^2))*lambda + sum(sum(B.^2))*lambda2 + lambda3*sum(sum(U.^2));
    

    Z = (Y*B' + lambda*F*U)*inv(B*B' + lambda*eye(size(B,1)));
    
    
    
    B = (Z'*Z + lambda2*eye(size(Z,2)))\Z'*Y;
    
    U = (F'*F*lambda + lambda3*eye(size(U,1)))\(lambda*F'*Z);
   
   if ii>1
      if (abs(o(ii)-o(ii-1))/o(ii))<tol
         break
      end
    end

end

erroro = sum(sum((Y-Z*B).^2))./sum(sum(Y.^2)) + sum(sum((Z-F*U).^2))./sum(sum((F*U).^2));
error1 = sum(sum((Y-Z*B).^2))./sum(sum(Y.^2)) ;
error2 = sum(sum((Z-F*U).^2))./sum(sum((F*U).^2));

dz = Z*(B*B' + lambda*eye(size(B,1)))-(Y*B' + lambda*F*U);
db = (Z'*Z + lambda2*eye(size(Z,2)))*B - Z'*Y;
du = (F'*F*lambda + lambda3*eye(size(U,1)))*U-lambda*F'*Z;





