function A = ml_csnnet_ssnet(X, idx)
% Cell-specific (single-sample) network degree matrix-based transformation
% https://academic.oup.com/nar/article/47/11/e62/5377474
% https://github.com/wys8c764/CSN
% https://github.com/WilfongGuo/Benchmark_control
    if nargin < 2, idx = 1; end
    pw1 = fileparts(mfilename('fullpath'));
    pth = fullfile(pw1, 'external', 'ml_CSN_transform');
    if ~(ismcc || isdeployed), addpath(pth); end
    A = csnet(X, idx);
end