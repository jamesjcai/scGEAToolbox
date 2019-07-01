function [X,Y]=run_hcp(X,Y)
% Run HCP (Hidden Covariates with Prior) to normalize RNA-seq data
% HCP (Mostafavi et al. 2013)
% doi: 10.1371/journal.pone.0068141

pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/HCP_hidden_covariates_with_prior');
addpath(pth);

n1=size(X,2);
n2=size(Y,2);
batch=[ones(1,n1) 2*ones(1,n2)]';
XY=[X Y]';

XYn = standardize(XY')';
Fn = standardize(batch')';

% (4) set the model parameters
k = 20;
lambda = 20;
lambda2 = 1;
lambda3 = 1;
iter = 100;

% (5) run HCP
[Z,B]=hidden_covariates_model(Fn,XYn,k,lambda,lambda2,lambda3,iter);

% (6) get the residual data:
Res=(XYn-Z*B)';

X=Res(:,1:n1);
Y=Res(:,n1+1:end);
