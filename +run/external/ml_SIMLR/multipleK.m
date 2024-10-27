function D_Kernels = multipleK(x)


N = size(x, 1);
KK = 0;
sigma = 2:-0.25:1;
Diff = (dist2(x));
[T, ~] = sort(Diff, 2);
[~, n] = size(Diff);
allk = 10:2:30;
t = 1;
for l = 1:length(allk)
    if allk(l) < (size(x, 1) - 1)
        TT = mean(T(:, 2:(allk(l) + 1)), 2) + eps;
        Sig = (repmat(TT, 1, n) + repmat(TT', n, 1)) / 2;
        Sig = Sig .* (Sig > eps) + eps;
        for j = 1:length(sigma)
            W = normpdf(Diff, 0, sigma(j)*Sig);
            Kernels(:, :, KK+t) = (W + W') / 2;
            t = t + 1;
        end
    end
end

for i = 1:size(Kernels, 3)
    K = Kernels(:, :, i);
    k = 1 ./ sqrt(diag(K)+1);
    %G = K.*(k*k');
    G = K;
    D_Kernels(:, :, i) = (repmat(diag(G), 1, length(G)) + repmat(diag(G)', length(G), 1) - 2 * G) / 2;
    D_Kernels(:, :, i) = D_Kernels(:, :, i) - diag(diag(D_Kernels(:, :, i)));

end

end