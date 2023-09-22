function [pc] = e_sparsepca(data, k, centerdata)

% data (rows: samples, columns: features)
if nargin < 3, centerdata = false; end
if centerdata
    mu = mean(data);
    data = bsxfun(@minus, data, mu);
end
[U, ~, ~] = svds(data', k);
pc = data * U;