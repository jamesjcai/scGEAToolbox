function [R] = e_distcorrmtx(X)

[n, m] = size(X);
m2 = m * m;

% Precompute centered distance matrices and variances
dvar = zeros(n, 1);
Avec = cell(n, 1);
for k = 1:n
    xk = X(k, :)';
    a = abs(xk - xk');
    a = a - mean(a) - mean(a, 2) + mean(a(:));
    dvar(k) = a(:)' * a(:) / m2;
    Avec{k} = a(:);
end

% Pairwise distance correlations — cache outer vector only
R = zeros(n);
for k = 1:n - 1
    ak = Avec{k};
    Avec{k} = [];  % free memory as we go
    for l = k + 1:n
        dcov = (ak' * Avec{l}) / m2;
        R(k, l) = sqrt(dcov / sqrt(dvar(k) * dvar(l)));
        R(l, k) = R(k, l);
    end
end
%{
any way to make this code faster:

function [R] = e_distcorrmtx(X)

n = size(X, 1); % number of genes
R = zeros(n);
c = 0;
for k = 1:n - 1
    fprintf('%d......%d\n', k, n)
    for l = k + 1:n
        c = c + 1;
        R(k, l) = pkg.e_distcorr(X(k, :)', X(l, :)');
        R(l, k) = R(k, l);
        % R(k,l)=pkg.e_bcdistcorr(X(k,:)',X(l,:)');
        % if c>20, break; end
    end
    % if c>20, break; end
end
%}