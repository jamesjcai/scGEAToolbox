function [a, fano, cv] = i_grpmean(x, c)

ugroups = unique(c);
n = numel(ugroups);
ng = size(x, 1);
a = zeros(ng, n);
fano = zeros(ng, n);
cv = zeros(ng, n);

needvar = nargout > 1;
needstd = nargout > 2;

for k = 1:n
    xk = x(:, c == ugroups(k));
    a(:, k) = full(mean(xk, 2));
    if needvar
        fano(:, k) = full(var(xk, 0, 2));
    end
    if needstd
        cv(:, k) = full(std(xk, 0, 2));
    end
end

if needvar
    fano = fano ./ max(a, eps);
end
if needstd
    cv = cv ./ max(a, eps);
end
end
