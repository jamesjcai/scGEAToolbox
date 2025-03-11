function [Y] = sc_ifft_dev(X)
% https://github.com/Sanofi-Public/PMCB-scGFT/blob/master/R/Core.R

    X = pkg.norm_libsize(X, 1e4);
    X = log1p(X);
    ft_mtx = fft(X, [], 2);
    
    [M, N] = size(ft_mtx);
    % ft_mtx stores FFT freq. components in rows
    % M = number of genes (also components)
    a = abs(ft_mtx);
    a = a./vecnorm(a, 2, 2);
    sigma = var(a, 1, 2);
    sigma(1) = 1;
    
    idx = randi([2, M], N, 1);
    A = normrnd(2, sqrt(sigma(idx)));
    
    oldft1 = ft_mtx(1,:);
    for k = 1:N
        ft_mtx(idx(k), k) = A(k)*ft_mtx(idx(k), k);
        if idx(k) ~= 1 && mod(M, 2)==0 && idx(k) ~= M/2 + 1
            ft_mtx(M - idx(k), k) = A(k)*ft_mtx(M - idx(k), k);
        end
    end
    assert(isequal(oldft1, ft_mtx(1,:)))

    Y = ifft(ft_mtx, [], 2);
    Y = real(Y);
    Y = exp(Y)-1;
    Y(Y<0) = 0;
    Y = round(Y);
    % Relative Percent Difference (RPD)
    % devs = mean(abs(Y - X));
end