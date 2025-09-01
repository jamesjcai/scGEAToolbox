%{
addpath('RT')
load example_data.mat

VIM = GENIE3(data);

% Indices of the genes that are used as candidate regulators
input_idx = [1, 2, 6, 10];
VIM2 = GENIE3(data, input_idx);


% Use Extra-Trees method
tree_method = 'ET';
% Number of randomly chosen candidate regulators at each node of a tree
K = 4;
% Number of trees per ensemble
ntrees = 50;
% Run the method with these settings
VIM3 = GENIE3(data, input_idx, tree_method, K, ntrees);

get_link_list(VIM)
%}