function A = xicornet(X, symmetric)
% Construct co-expression network using Chatterjee's xi correlation
%
% A = net.xicornet(X)
% A = net.xicornet(X, symmetric)
%
% X         - genes x cells expression matrix (log-normalised recommended)
% symmetric - if true, use symmetric xi (default: false)
% A         - genes x genes adjacency matrix of xi correlation values
%
% Chatterjee's xi detects nonlinear and non-monotonic dependence.
% ref: Chatterjee, JASA 2021. DOI:10.1080/01621459.2020.1758115

arguments
    X double
    symmetric(1,1) logical = false
end

X = double(full(X));
n = size(X, 1);
A = zeros(n);

if symmetric
    % Xi correlation is symmetric — compute upper triangle only
    parfor k = 1:n - 1
        row = zeros(1, n);
        for l = k + 1:n
            row(l) = pkg.e_xicor(X(k, :), X(l, :), true);
        end
        A(k, :) = row;
    end
    A = A + A';
else
    parfor k = 1:n
        row = zeros(1, n);
        for l = 1:n
            if k ~= l
                row(l) = pkg.e_xicor(X(k, :), X(l, :), false);
            end
        end
        A(k, :) = row;
    end
end
end
