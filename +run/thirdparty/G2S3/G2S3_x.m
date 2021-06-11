function [data_impute, W] = G2S3_x(data,varargin)

pw0=pwd;
pw1 = fileparts(mfilename('fullpath'));
cd(pw1);

addpath(fullfile('gaimc','graphs'));
addpath("gaimc");
run(fullfile('gspbox','gsp_start.m'));
run(fullfile('unlocbox','init_unlocbox.m'));

% Graph of smooth gene network signal inputed genes.
% 
% data format: cell(row) * genes(columns) array
%   'a' (default = 1)
%       bigger a -> bigger weights in gene networks
%   'b' (default = 1)
%       bigger b -> more dense W   
%   'mode'(default = 'logb')
%       'logb': with -log(degree_node) barrier to the learning of graph  
%        minimize_W sum(sum(W .* Z)) - a * sum(log(sum(W))) + b * ||W||_F^2/2 + c * ||W-W_0||_F^2/2
%
%
%       'l2p': add L2 penalty on the node degree
%        minimize_W sum(sum(W .* Z)) + a/2 * ||W||_F^2/2 
%                 + a/2 * ||sum(W)||_2^2 
%                 s.t. sum(sum(W)) == n
%   'tol'(default = 1e-4)
%       cut-off for neglectable weights   
%   'normalize' (default = 1)
%       whether or not normalize the cell library in the output data
%   'scale' (default = 1)
%       whether or not scale the gene expression level back
% output format: cell(row) * genes(columns) array

% set up default parameters   
a = 1; b = 1; type = 'logb';tol = 10^-4;self = 1; normalize = 1;scale = 1;
% get input parameters

for i=1:length(varargin)
    % mode
    if(strcmp(varargin{i},'mode'))
       type = lower(varargin{i+1});
    end
    % a 
    if(strcmp(varargin{i},'a'))
       a = lower(varargin{i+1});
    end
    % b
    if(strcmp(varargin{i},'b'))
       b = lower(varargin{i+1});
    end
    % tol
    if(strcmp(varargin{i},'epsi'))
       tol = lower(varargin{i+1});
    end
    if(strcmp(varargin{i}, 'normalize'))
      normalize = lower(varargin{i+1});
    end

    if(strcmp(varargin{i}, 'self'))
      self = lower(varargin{i+1});
    end

    if(strcmp(varargin{i}, 'scale'))
      normalize = lower(varargin{i+1});
    end
  
end

[~, nc] = size(data);

% log transformation


M = normc(data);

function W = affinity2graph(Z)
    switch type
        case 'logb'
            [W, stat] = gsp_learn_graph_log_degrees(Z, a, b);
        case 'l2p'
            [W, stat] = gsp_learn_graph_l2_degrees(Z, a);

    W(W<tol) = 0;
    end
end

% if nc < 1000
    Z = gsp_distanz(M).^2;
    W = affinity2graph(Z);
% else
    
% learning large graph - use flann
    
%      params_NN.use_flann = 1;
%  
%      sub_idx = randsample(nc,1000);
%      M_sub = M(:, sub_idx);
%      Z_sub = gsp_distanz(M_sub).^2;
%      W_sub = affinity2graph(Z_sub);
%      
%    
%     k = 10; % average edges per node
%     params_NN.use_flann = 1;
%     params_NN.k = 2 * k + 1;
%     [indx, indy, dist, ~, ~, ~, NN, Z_sorted] = gsp_nn_distanz(data, data, params_NN);
%     Z_sp = sparse(indx, indy, dist.^2, nc, nc, params_NN.k * nc * 2);
%     Z_sp = gsp_symmetrize(Z_sp, 'full');
%     Z_sorted = Z_sorted(:, 2:end).^2; 
%     theta = gsp_compute_graph_learning_theta(Z_sorted, k, 0, true);
%     params.edge_mask = Z_sp > 0;
%     params.fix_zeros = 1;
%     [W, ~] = gsp_learn_graph_log_degrees(Z_sp * theta, a, b, params);
    
% end 

switch self
  case 1
    W = bsxfun(@rdivide, W, sum(W,1));
    W = W + eye(nc);
    disp('adding self prob')
  case 0
    W = W;
    disp('only use neighbor')
end
  

W(W<0) = 0;

W = bsxfun(@rdivide, W, sum(W,1));
W = W^2;
data_impute = data * W;



if normalize ==1
  data_impute = bsxfun(@rdivide, data_impute, sum(data_impute,2))*mean(sum(data_impute,2));
  disp('normalized by library size')
else
  %data_impute = data_impute;
  disp('unormalized')
end

if scale ==1
  data_impute = data_impute./sum(data_impute);
  data_impute = data_impute.*sum(data);
end
cd(pw0);

end