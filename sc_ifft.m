function [Y] = sc_ifft(X)
% https://github.com/Sanofi-Public/PMCB-scGFT/blob/master/R/Core.R
% n = nsynth            Number of synthetic cells to be generated.
% k = ncpmnts           The number of components to modify during the
%                             synthesis process.

X = pkg.norm_libsize(X, 1e4);
X = log1p(X);
ft_mtx = fft(X, [], 2);

[M, N] = size(ft_mtx);
% ft_mtx stores FFT freq. components in rows
% M = number of genes (also? components)

a = abs(ft_mtx);
a = a./vecnorm(a, 2, 2);
sigma = var(a, 1, 2);
sigma(1) = 1;

k = randi([2, M], N, 1);
A = normrnd(2, sqrt(sigma(k)));

for cn = 1:N
    ft_mtx(k(cn), cn) = A(cn)*ft_mtx(k(cn), cn);
    %if k(cn) ~= 1 && mod(M, 2)==0 && k(cn) ~= M/2 + 1
        ft_mtx(M - k(cn), cn) = A(cn)*ft_mtx(M - k(cn), cn);
    %end
end

Y = ifft(ft_mtx, [], 2);

end
