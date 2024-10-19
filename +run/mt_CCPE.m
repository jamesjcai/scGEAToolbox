function [pseudotime] = mt_CCPE(X)

    pw1 = fileparts(mfilename('fullpath'));
    if ~(ismcc || isdeployed)
        addpath(fullfile(pw1, 'external', 'mt_CCPE'));
    end
    
    X = sc_transform(X,"type","PearsonResiduals");
    lambda=70;
    gamma=140;
    sigma=0.001;  %Gaussian distribution
    [pseudotime] = CCPE(X, lambda, gamma, sigma);

end