% Author: Van Hoan Do
function cls = eval_exact_Specter(fea, n_clusters, ensemble_size, mingamma)
% Input: 
    % -fea: expression data where rows are cells, collumns are principal components computed by PCA or genes
    % -n_clusters: number of clusters
    % -ensemble_size: number of clusterings in the ensemble 
    % -mingamma: minimum gaussion bandwidth (default: 0.1)
% Output:
    % cls: clusters of cells.
    
    % matlab seed
    rand ("state", 1000000);
    rand('seed', 1089); 
    [m, n] = size(fea);
    % apply pre-processing
    params.mode =0;
    params.HV = 1; % use highly variable genes selection: 1: yes, 0: no
    params.PCA = 1; % apply PCA: 1: yes, 0: no
    params.print = 1; % print results
    params.n_clusters = n_clusters;
    [opts, fea] = learn_LSC(fea, params);
    opts.seed = 0;
    
    N = ensemble_size;
    
    clusters = zeros(N,m);
    optsVec = repmat(opts, N);
    parfor i=1:1:N
        optsVec(i).r = max(5, round(opts.r*(rand*0.5+0.8)));
        optsVec(i).seed = 10*i; % use different seed to avoid repetive clustering.
        clusters(i,:) = LSC_eigen(fea, n_clusters, optsVec(i), rand*0.1 + mingamma); 
    end

    cls = evalCOAL(clusters, n_clusters); % run hierachical clustering of AL on co-association matrix

end
