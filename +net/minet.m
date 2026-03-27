function MI_mat = minet(X)
% Construct co-expression network using pairwise mutual information
%
% MI_mat = net.minet(X)
%
% X      - genes x cells expression matrix
% MI_mat - genes x genes mutual information matrix (symmetric)
%
% Uses parallel computation. Pairwise MI computed via binned histograms.
% See also: qtm.BinPairMI

X = double(full(X));

% Transpose for efficient column access across observations
data = X.';
ngene = size(data, 2);

MI_mat = zeros(ngene);
% Compute upper triangle in parallel; fill both sides directly
for jg = 1:ngene - 1
    tmp = zeros(ngene - jg, 1);
    d_jg = data(:, jg);
    parfor idx = 1:ngene - jg
        tmp(idx) = qtm.BinPairMI(d_jg, data(:, jg + idx));
    end
    MI_mat(jg, jg+1:ngene) = tmp';
    MI_mat(jg+1:ngene, jg) = tmp;
end
end
