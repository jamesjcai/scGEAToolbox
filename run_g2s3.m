function X=run_g2s3(X)
% G2S3: a gene graph-based imputation method

pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/G2S3');
addpath(pth);

% G2S3 needs [cells x genes]
X=X';
[X] = G2S3_x(X);
X=X';


