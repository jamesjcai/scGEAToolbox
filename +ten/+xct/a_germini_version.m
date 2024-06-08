% 
% 
% classdef ManifoldAlignmentNet
%     % This class defines a neural network for manifold alignment.
% 
%     properties
%         n_models;  % Number of data models
%         data_arr;  % Cell array containing data matrices
%         w;        % Sparse correspondence matrix
%         model_dic; % Dictionary containing trained models
%         projections; % Low-dimensional projections of data
%         losses;    % Training loss history
%     end
% 
%     methods
%         function obj = ManifoldAlignmentNet(data_arr, w, n_dim, layers)
%             % Constructor
%             [obj.n_models, obj.data_arr, obj.w] = obj._check_data(data_arr, w);
%             obj.model_dic = obj.create_models(layers, n_dim);
%         end
% 
%         function [n_models, data_arr, w] = _check_data(self, data_arr, w)
%             % Check data consistency
%             n_models = length(data_arr);
%             if ~all(isa(data_arr{:}, 'double'))
%                 error('Input a list of counts in double arrays with genes by cells');
%             end
%             if ~sum([size(x, 1) for x = data_arr]) == size(w, 1)
%                 error('Input sequence of counts consistent with correspondence');
%             end
%         end
% 
%         function model_dic = create_models(self, layers, n_dim)
%             % Create models with specified architecture
%             model_dic = containers.Map;
%             if isempty(layers)
%                 a = 4;
%                 for i = 1:self.n_models
%                     n_h = round(geomean([size(self.data_arr{i}, 2), n_dim]));
%                     model_dic(sprintf('model_%d', i)) = [a * n_h, n_h, n_dim];
%                 end
%             else
%                 for i = 1:self.n_models
%                     model_dic(sprintf('model_%d', i)) = layers;
%                 end
%             end
% 
%             % Initialize models and convert data to tensors
%             rng(0);  % Set random seed for reproducibility
%             for i = 1:self.n_models
%                 model_dic(sprintf('model_%d', i)) = Net(size(self.data_arr{i}, 2), model_dic(sprintf('model_%d', i)){:});
%                 self.data_arr{i} = gpuArray(single(self.data_arr{i}));  % Move data to GPU if available
%             end
%         end
% 
%         function arch(self)
%             % Print model architecture
%             for i = 1:self.n_models
%                 fprintf('model_%d:\n', i);
%                 disp(model_dic(sprintf('model_%d', i)));
%             end
%         end
% 
%         function save_model_states(self, file_dir)
%             % Save model states
%             mkdir(file_dir);
%             for i = 1:self.n_models
%                 save(fullfile(file_dir, sprintf('model_%d.mat', i)), 'model_dic', sprintf('model_%d', i));
%                 fprintf('Saved model to %s/model_%d.mat\n', file_dir, i);
%             end
%         end
% 
%         function load_model_states(self, file_dir)
%             % Load model states
%             for i = 1:self.n_models
%                 load(fullfile(file_dir, sprintf('model_%d.mat', i)), sprintf('model_%d', i));
%                 model_dic(sprintf('model_%d', i)).eval();
%                 fprintf('Loaded model from %s/model_%d.mat\n', file_dir, i);
%             end
%         end
% 
%         function projections = reload_embeds(self)
%             % Recompute embeddings
%             y_preds = cell(1, self.n_models);
%             for i = 1:self.n_models
%                 y_preds{i} = model_dic(sprintf('model_%d', i))(self.data_arr{i});
%             end
%             outputs = cat(1, y_preds{:});
%             [u, ~, v] = svd(outputs, 'econ');
%             proj_outputs = u * v';
