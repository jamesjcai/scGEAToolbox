% classdef a_chatgpt_version < torch.nn.Module
%     properties
%         linear1
%         linear2
%         linear3
%     end
%     methods
%         function obj = a_chatgpt_version(D_in, H1, H2, D_out)
%             obj.linear1 = torch.nn.Linear(D_in, H1);
%             obj.linear2 = torch.nn.Linear(H1, H2);
%             obj.linear3 = torch.nn.Linear(H2, D_out);
%         end
% 
%         function y_pred = forward(obj, x)
%             h1_sigmoid = obj.linear1(x).sigmoid();
%             h2_sigmoid = obj.linear2(h1_sigmoid).sigmoid();
%             y_pred = obj.linear3(h2_sigmoid);
%         end
%     end
% end
% 
% classdef ManifoldAlignmentNet
%     properties
%         n_models
%         data_arr
%         w
%         model_dic
%         projections
%         losses
%     end
% 
%     methods
%         function obj = ManifoldAlignmentNet(data_arr, w, n_dim, layers)
%             % Constructor
%             obj.n_models = length(data_arr);
%             obj.data_arr = data_arr;
%             obj.w = w;
%             obj.model_dic = obj.create_models(layers, n_dim);
%         end
% 
%         function model_dic = create_models(obj, layers, n_dim)
%             layer_dic = containers.Map;
%             if isempty(layers)
%                 a = 4;
%                 for i = 1:obj.n_models
%                     n_h = geomean([size(obj.data_arr{i}, 2), n_dim]);
%                     layer_dic(num2str(i)) = [a * n_h, n_h, n_dim];
%                 end
%             else
%                 for i = 1:obj.n_models
%                     layer_dic(num2str(i)) = layers;
%                 end
%             end
% 
%             model_dic = containers.Map;
%             rng(0); % Set random seed
%             for i = 1:obj.n_models
%                 model_dic(['model_', num2str(i)]) = Net(size(obj.data_arr{i}, 2), layer_dic(num2str(i))(1), layer_dic(num2str(i))(2), layer_dic(num2str(i))(3));
%                 obj.data_arr{i} = torch.tensor(obj.data_arr{i}, 'dtype', torch.float32);
%             end
%         end
% 
%         function projections = reload_embeds(obj)
%             y_preds = [];
%             for i = 1:obj.n_models
%                 y_preds = [y_preds; obj.model_dic(['model_', num2str(i)])(obj.data_arr{i})];
%             end
%             outputs = torch.cat(y_preds, 1);
%             [~, ~, v] = svd(outputs, 'econ');
%             proj_outputs = outputs * v';
%             obj.projections = proj_outputs.detach().numpy();
%             projections = obj.projections;
%         end
% 
%         function train(obj, n_steps, lr, verbose, optim_kwargs)
%             assert(n_steps > 0);
%             obj.losses = [];
%             L_np = full(csgraph.laplacian(obj.w, 'unnormalized'));
%             L = torch.tensor(L_np, 'dtype', torch.float32);
%             params = cell(1, obj.n_models);
%             for i = 1:obj.n_models
%                 params{i} = obj.model_dic(['model_', num2str(i)]).parameters();
%             end
%             optimizer = torch.optim.Adam([params{:}], 'lr', lr, optim_kwargs{:});
% 
%             for i = 1:obj.n_models
%                 obj.model_dic(['model_', num2str(i)]).train();
%             end
% 
%             for t = 1:n_steps
%                 y_preds = [];
%                 for i = 1:obj.n_models
%                     y_preds = [y_preds; obj.model_dic(['model_', num2str(i)])(obj.data_arr{i})];
%                 end
%                 outputs = torch.cat(y_preds, 1);
%                 [u, ~, v] = svd(outputs, 'econ');
%                 proj_outputs = u * v';
%                 loss = trace(proj_outputs' * L * proj_outputs) / 3000;
%                 obj.losses(end + 1) = loss.item();
% 
%                 if verbose || mod(t, 10) == 0
%                     disp([t, loss.item()]);
%                 end
% 
%                 optimizer.zero_grad();
%                 loss.backward();
%                 optimizer.step();
%             end
%             obj.projections = proj_outputs.detach().numpy();
%         end
% 
%         function dist_df = nn_aligned_dist(obj, projections, gene_names_x, gene_names_y, w12_shape, dist_metric, rank, verbose)
%             if verbose
%                 disp(['computing pair-wise ', dist_metric, ' distances...']);
%             end
%             dist = pdist2(projections(1:size(projections, 1) / 2, :), projections(size(projections, 1) / 2 + 1:end, :), dist_metric);
%             dist_df = array2table(dist, 'RowNames', gene_names_x, 'VariableNames', gene_names_y);
%             dist_df = stack(dist_df);
%             dist_df.Properties.VariableNames = {'ligand', 'receptor', 'dist'};
%             dist_df.Properties.RowNames = strcat(dist_df.ligand, '_', dist_df.receptor);
%             w12 = obj.w(1:w12_shape(1), w12_shape(2):end);
%             dist_df.correspondence = reshape(w12, [numel(w12), 1]);
%             if rank
%                 dist_df = sortrows(dist_df, 'dist');
%                 dist_df.rank = (1:size(dist_df, 1))';
%             end
%         end
% 
%         function dist_df_filtered = filtered_nn_aligned_dist(~, df_nn, candidates)
%             if ismember('diff2', df_nn.Properties.VariableNames)
%                 df_nn_filtered = df_nn(candidates, :);
%                 df_nn_filtered = sortrows(df_nn_filtered, 'diff2');
%             else
%                 df_nn_filtered = df_nn(candidates, :);
%                 df_nn_filtered = sortrows(df_nn_filtered, 'dist');
%             end
%             disp(['manifold aligned # of L-R pairs:', num2str(size(df_nn_filtered, 1))]);
%             df_nn_filtered.rank_filtered = (1:size(df_nn_filtered, 1))';
%             dist_df_filtered = df_nn_filtered;
%         end
%     end
% end
% 
