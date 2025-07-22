function [c] = ml_Specter(X, k)

if nargin < 2, k = 5; end

if ~(ismcc || isdeployed)
    pw1 = fileparts(mfilename('fullpath'));
    pth0 = fullfile(pw1, '..', 'external', 'ml_Specter');
    pth1 = fullfile(pw1, '..', 'external', 'ml_Specter', 'dimred');
    pth2 = fullfile(pw1, '..', 'external', 'ml_Specter', 'LSC');
    pth3 = fullfile(pw1, '..', 'external', 'ml_Specter', 'utils');
    addpath(pth0);
    addpath(pth1);
    addpath(pth2);
    addpath(pth3);
end
data = transpose(sc_transform(X));

%% choose parameters for algorithm
n_clusters = k; % number of clusters
ensemble_size = 200; % ensemble sizes (default: 200)
mingamma = 0.1; % minimum gaussian bandwidth (default: 0.1), best value in range 0.1 to 0.8
% All gamma value will be chosen from the inteval [mingamma, mingamma + 0.1]
% Please refer to the paper for the details

%% We have developed two version of Specter. Let Specter to decide the algorithm
c = eval_auto_Specter(data, n_clusters, ensemble_size, mingamma);

end
