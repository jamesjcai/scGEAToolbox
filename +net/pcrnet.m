function [A] = pcrnet(X, ncom, fastersvd, dozscore, UseParallel, guiwaitbar)
% Construct GRN using principal component regression (PCR)
%
% A = net.pcrnet(X)
% A = net.pcrnet(X, ncom)
% A = net.pcrnet(X, ncom, fastersvd, dozscore, UseParallel, guiwaitbar)
%
% X           - genes x cells expression matrix (LogNormalized recommended)
% ncom        - number of principal components (default: 3)
% fastersvd   - use lmsvd instead of svds (default: true)
% dozscore    - z-score genes before regression (default: true)
% UseParallel - use parfor loop (default: true)
% guiwaitbar  - show GUI progress bar, serial only (default: false)
%
% ref: https://rdrr.io/cran/dna/man/PCnet.html
%      https://github.com/cran/dna/blob/master/src/rpcnet.c

arguments
    X double
    ncom(1, 1) {mustBeNumeric} = 3
    fastersvd(1, 1) logical = false
    dozscore(1, 1) logical = true
    UseParallel(1, 1) logical = false
    guiwaitbar(1, 1) logical = false
end

opts.maxit = 150;

% Ensure toolbox root is on worker path so +net package is visible
if UseParallel
    rootdir = fileparts(fileparts(mfilename('fullpath')));
    addpath(rootdir);
end

X = X.';
if dozscore
    X = zscore(X);
end
n = size(X, 2);
A = 1 - eye(n);

% Precompute column index mask to avoid copying X and deleting a column
idx_all = 1:n;

if UseParallel
    B = A(:, 1:end-1);
    warning off
    parfor k = 1:n
        y = X(:, k);
        cols = [idx_all(1:k-1), idx_all(k+1:end)];
        Xi = X(:, cols);
        if fastersvd
            [~, ~, coeff] = lmsvd(Xi, ncom, opts);
        else
            [~, ~, coeff] = svds(Xi, ncom);
        end
        score = Xi * coeff;
        nrm2 = sum(score .* score);
        Beta = sum(y .* (score ./ nrm2));
        B(k, :) = coeff * Beta';
    end
    warning on
    for k = 1:n
        A(k, A(k, :) == 1) = B(k, :);
    end
else
    if guiwaitbar
        fw = gui.gui_waitbar_adv;
    end
    for k = 1:n
        if guiwaitbar
            gui.gui_waitbar_adv(fw, k/n);
        end
        y = X(:, k);
        cols = [idx_all(1:k-1), idx_all(k+1:end)];
        Xi = X(:, cols);
        if fastersvd
            [~, ~, coeff] = lmsvd(Xi, ncom, opts);
        else
            [~, ~, coeff] = svds(Xi, ncom);
        end
        score = Xi * coeff;
        nrm2 = sum(score .* score);
        Beta = sum(y .* (score ./ nrm2));
        A(k, A(k, :) == 1) = coeff * Beta';
    end
    if guiwaitbar, gui.gui_waitbar_adv(fw); end
end
end
