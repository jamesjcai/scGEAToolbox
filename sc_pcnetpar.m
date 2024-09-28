function [A] = sc_pcnetpar(X, ncom, fastersvd, dozscore)
% [A]=sc_pcnetpar(X,ncom)
% ncom - number of components used (default=3)
% ref: https://rdrr.io/cran/dna/man/PCnet.html
% https://github.com/cran/dna/blob/master/src/rpcnet.c
% https://rdrr.io/cran/dna/f/inst/doc/Introduction.pdf

if nargin < 4, dozscore = true; end
if nargin < 3, fastersvd = true; end
if nargin < 2, ncom = 3; end

opts.maxit = 150;

% if fastersvd
%     pw1 = fileparts(mfilename('fullpath'));
%     pth = fullfile(pw1, '+run', 'thirdparty', 'faster_svd', 'lmsvd');
%     if ~(ismcc || isdeployed), addpath(pth); end
% end

% X is supposed to be LogNormalized, i.e., [X]=log1p(sc_norm(X));
X = X';
if dozscore
    X = zscore(X);
end

n = size(X, 2);
A = 1 - eye(n);
B = A(:, 1:end-1);
warning off
parfor k = 1:n
    y = X(:, k);
    Xi = X;
    Xi(:, k) = [];
    if fastersvd
        % disp('Using fastsvd.')
        [~, ~, coeff] = lmsvd(Xi, ncom, opts);
    else
        [~, ~, coeff] = svds(Xi, ncom);
        %[~,~,coeff]=rsvd(Xi,ncom);
    end
    score = Xi * coeff;
    score = (score ./ (vecnorm(score).^2));
    Beta = sum(y.*score);
    B(k, :) = coeff * Beta';
end
warning on
for k = 1:n
    A(k, A(k, :) == 1) = B(k, :);
end
end