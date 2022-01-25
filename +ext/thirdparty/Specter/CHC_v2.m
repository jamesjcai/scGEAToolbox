% function [C_est, lk_est, time, IDX_LD, ind_obs, weight_VD] = CSC(G,param_CSC)
%
% This is the main function of the CSC toolbox. 
% 
% Given a graph structure G where:
% - G.W contains the sparse weighted symmetrical adjacency matrix of the graph
% - G.k is the number of classes you are looking for
% - G.N is the number of nodes of the graph (should be equal to length(G.W)
% 
% Given parameters param_CSC where: 
% - param_CSC.poly_order is the order of the Jackson-chebychev polynomial approximation. 
% Default is 50.
% - param_CSC.regu is the regularisation parameter of the decoder used for interpolation. 
% Default is 1e-3.
% - param_CSC.sampling is the sampling distribution you want. Either 'uniform' or 'VD' 
% for variable density sampling. Default is 'uniform'. 
% - param_CSC.n_factor is the factor for the sampled number of nodes. The number of sampled nodes 
% reads: n = param_CSC.n_factor x G.k * log(G.k). Default is 2.
% - param_CSC.d_factor is the factor for the number of random vectors. The number of random vectors 
% reads: d = param_CSC.d_factor x log(n). Default is 4.
% - param_CSC.lap_type is the choosen Laplacian. Either 'combinatorial' or 'normalized'. 
% Default is 'normalized'.
% - param_CSC.solver is the MATLAB solver used for interpolation. Either 'gmres' or 'cgs'. 
% Default is 'gmres'
% If some (or all) fields of param_CSC are not given, the function uses the default values.  
% 
% This function outputs:
% - C_est is the G.N x G.k matrix, where column j is the recovered indicator vector of class j. 
% It still needs to be normalized and binarized to obtain a strict partition of the graph. 
% - lk_est the estimated k-th eigenvalue of its Laplacian
% - time is a structure where are recorded the computation times of the different steps of the algorithm
% - IDX_LD is the result of the low-dimensional k-means
% - ind_obs are the indices of the n sampled nodes
% - weight_VD is the estimated probability distribution if one wants to use variable density sampling
function [C_est, lk_est, time, IDX_LD, ind_obs, weight_VD] = CSC_v2(G,param_CSC,IDX_LD, ind_obs)
%%%
%% ====================== check parameter list ================================
%%%

if ~isfield(param_CSC, 'poly_order'),
    param_CSC.poly_order = 50;
end

if ~isfield(param_CSC, 'regu'),
    param_CSC.regu = 1e-3;
end

if ~isfield(param_CSC, 'sampling'),
    param_CSC.sampling = 'uniform';
end

if ~isfield(param_CSC, 'n_factor'),
    % n = param_CSC.n_factor * G.k * log(G.k);
    param_CSC.n_factor = 2;
end

if ~isfield(param_CSC, 'd_factor'),
    % d = param_CSC.d_factor * log(n);
    param_CSC.d_factor = 4;
end

if ~isfield(param_CSC, 'lap_type'),
    param_CSC.lap_type = 'normalized';
end

if ~isfield(param_CSC, 'solver'),
    param_CSC.solver = 'gmres';
end

%%%
%% ==================== compute required Laplacian matrix ===========================
%%%
if strcmp(param_CSC.lap_type,'normalized')
    normBOOL=true;
    G.lap_type='normalized';
    G.L = param_CSC.lap;
elseif strcmp(param_CSC.lap_type,'combinatorial')
    normBOOL=false;
    G.lap_type='combinatorial';
    G.L = param_CSC.lap;
    % G.L = create_laplacian(G.W,G.lap_type);
else 
    error(['CSC.m : param_CSC.lap_type should be ''normalized'' or ''combinatorial''']);
end

%%%
%% ====================== estimate lambda_k and weight VD ===========================
%%%
% normBOOL = 1
if normBOOL
    G.lmax=2;
    time.lmax=0;
else 
    tic;
    opts.isreal=1;opts.issym=1;opts.maxit=1000; %10000;
    G.lmax=eigs(G.L,1,'LA',opts);
    toc;
    % time.lmax=toc;
end

fprintf('Estimating lambda_k...')
tic;
param.order=param_CSC.poly_order;
if strcmp(param_CSC.lap_type,'normalized')
    [~, lk_est, cum_coh_k, ~] = estimate_lambda_k(G, G.k, param);
    % lk_est = 1.5 
else
    [~, lk_est, cum_coh_k, ~] = estimate_lambda_k(G, G.k, param);
    if normBOOL 
        lk_est = 1.0 
    end
end 
% 
% param.hint_lambda_max=lk_est*2;   
% [~, lk_estp1, cum_coh_k_p1, ~] = estimate_lambda_k(G, G.k, param);
% 
% lk_est=(lk_est+lk_estp1)/2;
% time.lk_est=toc;
% mean_num_coh=mean([cum_coh_k,cum_coh_k_p1],2);
% weight_VD = sum(mean_num_coh,2); 
% weight_VD=weight_VD./sum(weight_VD); 
weight_VD = 0;
toc;

G.lk=(G.lmax + lk_est)/2;
n=length(ind_obs);
if strcmp(param_CSC.sampling,'uniform')
    weight = ones(G.N, 1)/(G.N); % Uniform density
elseif strcmp(param_CSC.sampling,'VD')
    weight = weight_VD; % Variable density
else error(['CSC: param_CSC.sampling must be either set to ''uniform'' or to ''VD''']);
end

% ind_obs = datasample(1:G.N, n, 'Replace', false, 'Weights', weight);

%%%
%% ====================== do k-means in low dimension ===========================
%%%

% IDX_LD = litekmeans(X_lk_est_DS , G.k,'MaxIter',1000, 'Replicates', 20);
% IDX_LD = runLWEA(G.W(ind_obs, ind_obs), G.k);
% time.k_means_low_dim=toc;

%%%
%% ====================== Interpolate in high dimensions: ===========================
%%%

fprintf('Interpolation of cluster indicators...\n\n')
tic;
C_obs_LD = sparse(1:n, IDX_LD, 1, n, G.k);
[~,JCH_HP] = jackson_cheby_poly_coefficients(G.lk,G.lmax,[0,G.lmax],param_CSC.poly_order);

C_est = zeros(G.N, size(C_obs_LD,2));
if (G.N < 10000)
    parfor k=1:size(C_obs_LD,2)
       c_obs = C_obs_LD(:,k);
       C_est(:, k) = interpolate_on_complete_graph(c_obs, ind_obs, @(x)gsp_cheby_op(G, JCH_HP, x), param_CSC.regu, G.N, param_CSC.solver);
    end
else 
    for k=1:size(C_obs_LD,2)
       c_obs = C_obs_LD(:,k);
       C_est(:, k) = interpolate_on_complete_graph(c_obs, ind_obs, @(x)gsp_cheby_op(G, JCH_HP, x), param_CSC.regu, G.N, param_CSC.solver);
    end
end
toc;
% time.interpolation=toc;
% time.total=time.interpolation+time.lmax+time.lk_est+time.k_means_low_dim;
time.total = 0;
