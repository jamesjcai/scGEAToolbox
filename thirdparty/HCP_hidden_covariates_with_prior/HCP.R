hidden_covariate_linear <- function(F,Y,k,lambda,lambda2,lambda3,iter=100,tol = 1e-6){
###########################################################################
# [Z,B,U,o,error,error1,error2,dz,db,du] = hidden_covariate_linear(F,Y,k,lambda,lambda2,lambda3,iter)
# Inputs:
#     F: a matrix nxd of known covariates, where n is the number of subjects and d is the number of known covariates. * Must not be mean centered/standardized
#     Y: a matrix of nxg of expression data (must be standardized (columns scaled to have constant SS and mean 0). * Must not be mean centered/standardized
#     k: number of inferred hidden components (k is an integer) lambda, lambda2, lambda3 are model parameters
#     iter (optional): number of iterations (default = 100)
#     tol (optional): convergence tolerance
#     note: k>0, lambda>0, lambda2>0, lambda3>0 must be set by the user based on the data at hand. one can set these values using cross-validation, by evaluating the "performance" of the  resulting residual data on a desired task.
#     typically, if lambda>5, then hidden factors match the known covariates closely.
#
# Objective:
#
#     this function solves the following problem:
#           argmin_{Z,B,U}   ||Y-Z*B||_2 + \lambda*||Z-F*U||_2 + \lambda2*||B||_2 + \lambda_3||U||_2
#
# Outputs: list consisting of following elements
#      R: matrix of residuals corrected for known and hidden covariates, dimesionality: nxg
#      Z: matrix of hidden components, dimensionality: nxk
#      B: matrix of effects of hidden components, dimensionality: kxg
#      o: value of objective function on consecutive iterations.
#
# to use the residual data: Residual = Y - Z*B
#
############################################################################################################

#### Standardize/Mean center inputs ####
Yn = scale(Y)
Fn = scale(F)

#### Set default values for outputs and initialise B ####
U = matrix(0,dim(F)[2],k)
Z = matrix(0,dim(F)[1],k)
B = matrix(runif(k*dim(Y)[2]),k,dim(Y)[2])

#### Error checking ####
if (dim(F)[1] != dim(Y)[1])
  stop('number of rows in F and Y must agree')

if (k<1 | lambda<tol | lambda2<tol | lambda3<tol)
  stop('lambda, lambda2, lambda3 must be positive and/or k must be an integer');

#### Coefficients esitmation ####
o = matrix(0,1,iter)
for (ii in 1:iter){
  o[ii] = norm(Yn-Z %*% B,type='E') + norm(Z-Fn %*% U, type='E')*lambda + norm(B,type='E')*lambda2 + norm(U,type='E')*lambda3
  
  Z = ((Yn %*% t(B)) + lambda * (Fn %*% U)) %*% solve((B %*% t(B)) + lambda*diag(dim(B)[1]))
  
  B = solve(((t(Z) %*% Z) + lambda2 * diag(dim(Z)[2])),(t(Z) %*% Yn))
  
  U = solve(((t(Fn) %*% Fn) * lambda + lambda3 * diag(dim(U)[1])),(lambda* t(Fn) %*% Z))
  
  if (ii > 1)
    if ((abs(o[ii]-o[ii-1])/o[ii]) < tol)
      break 
}

if (ii>=iter)
  writeLines("\nPre-mature convergence: Consider increasing the number of iterations")

#### Error calculation ####
error = (norm((Yn-Z %*% B),type='E') /norm(Yn,type='E')) + (norm((Z - Fn %*% U),type='E')/norm(Fn%*%U,type='E'))
error1 = norm(Yn-Z%*%B,type='E')/norm(Yn,type='E')
error2 = norm(Z-Fn%*%U,type='E')/norm(Fn%*%U,type='E')

#### Delta change calculation ####
dz = Z %*% (B %*% t(B) + lambda * diag(dim(B)[1])) - (Yn %*% t(B) + lambda * Fn %*% U)
db = (t(Z) %*% Z + lambda2 * diag(dim(Z)[2])) %*% B - t(Z) %*% Yn
du = (t(Fn) %*% Fn *lambda + lambda3 * diag(dim(U)[1])) %*% U - lambda * (t(Fn) %*% Z)

#### Residual calculation and covariate adjustment ####
R = (Yn - Z %*% B) * t(replicate(dim(Y)[1],attr(Yn,'scaled:scale'))) + t(replicate(dim(Y)[1],attr(Yn,'scaled:center')))

#### Final Outputs ####
return(list(R=R,Z=Z,B=B,U=U,o=o,error=error,error1=error1,error2=error2,dz=dz,db=db,du=du))
}