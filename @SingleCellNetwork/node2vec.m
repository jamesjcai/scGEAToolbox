function [embedding] = node2vec(A, d, walks_per_node, walk_length, p, q)
% NODE2VEC computes node embeddings using the Node2Vec algorithm
%
% Input:
%   A - adjacency matrix of the graph
%   d - dimensionality of the embedding
%   walks_per_node - number of walks to generate per node
%   walk_length - length of each random walk
%   p - parameter controlling likelihood of returning to previous node
%   q - parameter controlling likelihood of exploring new nodes
%
% Output:
%   embedding - n x d matrix, where n is the number of nodes in the graph
%               and d is the dimensionality of the embedding

n = size(A, 1); % number of nodes in the graph

% precompute transition probabilities
P = compute_transition_probs(A, p, q);

% generate random walks
walks = generate_walks(A, walks_per_node, walk_length, P);

% learn embeddings using Skip-Gram model
model = fit_skip_gram(walks, n, d);

% extract embeddings
embedding = model.W';
end

function [P] = compute_transition_probs(A, p, q)
% COMPUTE_TRANSITION_PROBS computes transition probabilities for each node
%
% Input:
%   A - adjacency matrix of the graph
%   p - parameter controlling likelihood of returning to previous node
%   q - parameter controlling likelihood of exploring new nodes
%
% Output:
%   P - transition probability matrix

n = size(A, 1); % number of nodes in the graph

% compute transition probabilities for each node
P = zeros(n);
for i = 1:n
    % get neighbors of current node
    nbrs = find(A(i, :));
    % compute unnormalized transition probabilities
    p_i = zeros(1, n);
    p_i(nbrs) = 1 ./ length(nbrs);
    for j = nbrs
        % compute second-order transition probabilities
        nbrs_j = find(A(j, :));
        p_ij = zeros(1, n);
        p_ij(nbrs_j) = q ./ length(nbrs_j);
        p_ij(i) = p;
        % combine first- and second-order probabilities
        p_i = p_i + p_ij;
    end
    % normalize probabilities
    P(i, :) = p_i ./ sum(p_i);
end
end

function [walks] = generate_walks(A, walks_per_node, walk_length, P)
% GENERATE_WALKS generates random walks on the graph
%
% Input:
%   A - adjacency matrix of the graph
%   walks_per_node - number of walks to generate per node
%   walk_length - length of each random walk
%   P - transition probability matrix
%
% Output:
%   walks - cell array of random walks

n = size(A, 1); % number of nodes in the graph

% generate random walks for each node
walks = cell(n*walks_per_node, 1);
for i = 1:n
    for j = 1:walks_per_node
        walk = zeros(1, walk_length);
        walk(1) = i;
        for k = 2:walk_length
            % sample next node from transition probabilities
            p_i = P(walk(k-1), :);
            walk(k) = randsample(n, 1, true, p_i);
        end
        walks{(i - 1)*walks_per_node+j} = walk;
    end
end
end


function [model] = fit_skip_gram(walks, n, d)
% FIT_SKIP_GRAM learns node embeddings using the Skip-Gram model
%
% Input:
%   walks - cell array of random walks
%   n - number of nodes in the graph
%   d - dimensionality of the embedding
%
% Output:
%   model - structure containing learned embeddings and other information

% concatenate all random walks into a single list of nodes
all_nodes = [walks{:}];

% count occurrences of each node
counts = histc(all_nodes, 1:n);

% compute unigram distribution
P = counts ./ sum(counts);

% initialize embeddings
W = randn(d, n);

% learn embeddings using Skip-Gram model
for i = 1:numel(walks)
    walk = walks{i};
    for j = 1:length(walk)
        % randomly subsample nodes to speed up training
        if rand < 1 - sqrt(0.001./P(walk(j)))
            continue;
        end
        % randomly select context nodes
        context = walk(max(1, j-5):min(length(walk), j+5));
        for k = 1:length(context)
            if context(k) == walk(j)
                continue;
            end
            % compute score for context node
            score = W(:, context(k))' * W(:, walk(j));
            % compute error
            err = 1 ./ (1 + exp(-score)) - (k == 1);
            % update embeddings
            W(:, walk(j)) = W(:, walk(j)) + 0.01 * err * W(:, context(k));
            W(:, context(k)) = W(:, context(k)) + 0.01 * err * W(:, walk(j));
        end
    end
end

model.W = W;
model.P = P;
end
