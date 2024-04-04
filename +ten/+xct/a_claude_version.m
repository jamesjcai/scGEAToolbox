Here's a translation of the provided Python code to MATLAB:

```matlab
% Import required libraries
import itertools
import numpy as np
import scipy
from scipy.sparse import coo_matrix
import torch
import torch.nn as nn

import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

from scTenifoldXct.stiefel import *
cuda = torch.device('cuda');

% Define the neural network
classdef Net < matlab.System & nnet.layer.Layer
    
    properties
        linear1
        linear2
        linear3
    end
    
    methods
        function layer = Net(D_in, H1, H2, D_out)
            layer.linear1 = nnet.layer.FullyConnectedLayer(H1, 'Activation', 'sigmoid');
            layer.linear2 = nnet.layer.FullyConnectedLayer(H2, 'Activation', 'sigmoid');
            layer.linear3 = nnet.layer.FullyConnectedLayer(D_out);
        end
        
        function y_pred = predict(layer, x)
            h1_sigmoid = layer.linear1.predict(x);
            h2_sigmoid = layer.linear2.predict(h1_sigmoid);
            y_pred = layer.linear3.predict(h2_sigmoid);
        end
    end
end

% Trainer class
classdef ManifoldAlignmentNet
    
    properties
        n_models
        data_arr
        w
        model_dic
    end
    
    methods
        function obj = ManifoldAlignmentNet(data_arr, w, n_dim, layers)
            [obj.n_models, obj.data_arr, obj.w] = obj.checkData(data_arr, w);
            obj.model_dic = obj.createModels(layers, n_dim);
        end
        
        function [n_models, data_arr, w] = checkData(obj, data_arr, w)
            n_models = length(data_arr);
            if ~all(cellfun(@(x) isa(x, 'double'), data_arr))
                error('Input a list of counts in numpy arrays with genes by cells');
            end
            if sum(cellfun(@(x) size(x, 1), data_arr)) ~= size(w, 1)
                error('Input sequence of counts consistent with correspondence');
            end
        end
        
        function model_dic = createModels(obj, layers, n_dim)
            layer_dic = containers.Map();
            if isempty(layers)
                a = 4;
                for i = 1:obj.n_models
                    n_h = floor(gmean([obj.data_arr{i}, n_dim]));
                    layer_dic(i) = {a * n_h, n_h, n_dim};
                end
            else
                for i = 1:obj.n_models
                    layer_dic(i) = layers;
                end
            end
            
            model_dic = containers.Map();
            rng(0); % Set seed
            for i = 1:obj.n_models
                model_dic(sprintf('model_%d', i)) = Net(size(obj.data_arr{i}, 2), layer_dic(i){:});
                obj.data_arr{i} = single(obj.data_arr{i});
            end
        end
        
        function arch(obj)
            for i = 1:obj.n_models
                fprintf('model_%d:\n', i);
                disp(obj.model_dic(sprintf('model_%d', i)));
            end
        end
        
        function saveModelStates(obj, file_dir)
            if ~exist(file_dir, 'dir')
                mkdir(file_dir);
            end
            for i = 1:obj.n_models
                save(sprintf('%s/model_%d.mat', file_dir, i), 'net');
                fprintf('save model to %s/model_%d.mat\n', file_dir, i);
            end
        end
        
        function loadModelStates(obj, file_dir)
            for i = 1:obj.n_models
                load(sprintf('%s/model_%d.mat', file_dir, i), 'net');
                obj.model_dic(sprintf('model_%d', i)) = net;
                fprintf('load model from %s/model_%d.mat\n', file_dir, i);
            end
        end
        
        function projections = reloadEmbeds(obj)
            y_preds = cell(1, obj.n_models);
            for i = 1:obj.n_models
                y_preds{i} = obj.model_dic(sprintf('model_%d', i)).predict(obj.data_arr{i});
            end
            outputs = vertcat(y_preds{:});
            [u, ~, v] = svd(outputs, 'econ');
            proj_outputs = u * v';
            obj.projections = double(proj_outputs);
        end
        
        function projections = train(obj, n_steps, lr, verbose, optim_kwargs)
            if nargin < 2
                n_steps = 1000;
            end
            if nargin < 3
                lr = 0.01;
            end
            if nargin < 4
                verbose = false;
            end
            if nargin < 5
                optim_kwargs = struct();
            end
            
            assert(n_steps > 0);
            obj.losses = [];
            L_np = scipy.sparse.csgraph.laplacian(obj.w, 'normed', false);
            L = single(L_np);
            
            params = cellfun(@(x) x.Learnables.Value, values(obj.model_dic), 'UniformOutput', false);
            optimizer = adam(params, 'MaxEpochs', n_steps, 'InitialLearnRate', lr, optim_kwargs{:});
            
            for i = 1:obj.n_models
                obj.model_dic(sprintf('model_%d', i)).Mode = 'train';
            end
            
            for t = 1:n_steps
                y_preds = cell(1, obj.n_models);
                for i = 1:obj.n_models
                    y_preds{i} = obj.model_dic(sprintf('model_%d', i)).predict(obj.data_arr{i});
                end
                
                outputs = vertcat(y_preds{:});
                [u, ~, v] = svd(outputs, 'econ');
                proj_outputs = u * v';
                
                loss = trace(proj_outputs' * L * proj_outputs) / 3000;
                
                if t == 1 || mod(t, 10) == 0
                    if verbose
                        fprintf('%d %f\n', t, loss);
                    end
                end
                obj.losses = [obj.losses, loss]; %#ok<AGROW>
                
                proj_outputs.LearnRateFactor = 1;
                
                [obj.model_dic.Values.LearnRateFactor] = deal(0);
                
                for i = 1:obj.n_models
                    model = obj.model_dic(sprintf('model_%d', i));
                    model.Learnables.LearnRateFactor = 1;
                end
                
                rgrad = proj_stiefel(proj_outputs, proj_outputs.Gradient);
                
                for i = 1:obj.n_models
                    model = obj.model_dic(sprintf('model_%d', i));
                    model.Learnables.Gradient = rgrad;
                end
                
                [obj.model_dic.Values.LearnRateFactor] = deal(0);
                
                [state, optim_info] = update(optimizer);
                
                for i = 1:obj.n_models
                    obj.model_dic(sprintf('model_%d', i)).Learnables.Value = state(i,:)';
                end
            end
            
            obj.projections = double(proj_outputs);
        end
        
        function plotLosses(obj, file