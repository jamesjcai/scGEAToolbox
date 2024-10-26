function X = ml_csndm_trans(X)
% Cell-specific (single-sample) network degree matrix-based transformation
% https://academic.oup.com/nar/article/47/11/e62/5377474
% https://github.com/wys8c764/CSN
% https://github.com/WilfongGuo/Benchmark_control
    pw1 = fileparts(mfilename('fullpath'));
    pth = fullfile(pw1, 'external', 'ml_CSN_transform');
    if ~(ismcc || isdeployed), addpath(pth); end
    % G=csnet(X(:,1:10));
    % csnedge(X(1,:),X(2,:));
    X = csndm(X);
end