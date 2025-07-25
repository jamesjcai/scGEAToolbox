function new_embedding = init_transform(indices, weights, embedding)
%INIT_TRANSFORM Given indices and weights and an original embedding
% initialize the positions of new points relative to the
% indices and weights (of their neighbors in the source data).
%
% new_embedding = INIT_TRANSFORM(indices, weights, embedding)
%
% Parameters
% ----------
% indices: array of size (n_new_samples, n_neighbors)
%     The indices of the neighbors of each new sample
% 
% weights: array of size (n_new_samples, n_neighbors)
%     The membership strengths of associated 1-simplices
%     for each of the new samples.
% 
% embedding: array of size (n_samples, dim)
%     The original embedding of the source data.
% 
% Returns
% -------
% new_embedding: array of size (n_new_samples, dim)
%     An initial embedding of the new sample points.

%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause
    
    new_embedding = zeros(size(indices,1), size(embedding,2));
    for i = 1:size(indices,1)
        for j = 1:size(embedding,2)
            for k = 1:size(indices,2)
                new_embedding(i, j) = new_embedding(i, j) + weights(i,k)*embedding(indices(i,k), j);
            end
        end
    end
end