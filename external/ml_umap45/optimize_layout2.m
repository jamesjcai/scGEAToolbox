function embedding = optimize_layout2(head_embedding, tail_embedding, ...
    head, tail, n_epochs, n_vertices, epochs_per_sample, a, b,...
    gamma, initial_alpha, negative_sample_rate, verbose, ~)
%OPTIMIZE_LAYOUT2 Another version of optimize_layout.m that slightly
% optimizes the speed for MATLAB; in particular, several quantities are
% vectorized. This function is only called if the UMAP method is 'MATLAB
% vectorized'.
%
% embedding = OPTIMIZE_LAYOUT2(head_embedding, tail_embedding, head, tail,
% n_epochs, n_vertices, epochs_per_sample, a, b)
%
% Parameters
% ----------
% head_embedding: array of size (n_samples, n_components)
%     The initial embedding to be improved by SGD.
% 
% tail_embedding: array of size (source_samples, n_components)
%     The reference embedding of embedded points. If not embedding new
%     previously unseen points with respect to an existing embedding this
%     is simply the head_embedding (again); otherwise it provides the
%     existing embedding to embed with respect to.
% 
% head: array of size (n_1_simplices, 1)
%     The indices of the heads of 1-simplices with non-zero membership.
% 
% tail: array of size (n_1_simplices, 1)
%     The indices of the tails of 1-simplices with non-zero membership.
% 
% n_epochs: double
%     The number of training epochs to use in optimization.
% 
% n_vertices: double
%     The number of vertices (0-simplices) in the dataset.
% 
% epochs_per_samples: array of size (n_1_simplices, 1)
%     A double value of the number of epochs per 1-simplex. 1-simplices with
%     weaker membership strength will have more epochs between being sampled.
% 
% a: double
%     Parameter of differentiable approximation of right adjoint functor.
% 
% b: double
%     Parameter of differentiable approximation of right adjoint functor.
% 
% gamma: double (optional, default 1)
%     Weight to apply to negative samples.
% 
% initial_alpha: double (optional, default 1)
%     Initial learning rate for the SGD.
% 
% negative_sample_rate: double (optional, default 5)
%     Number of negative samples to use per positive sample.
% 
% verbose: boolean (optional, default false)
%     Whether to report information on the current progress of the algorithm.
% 
% Returns
% -------
% embedding: array of size (n_samples, n_components)
%     The optimized embedding.
%
% See also: OPTIMIZE_LAYOUT

%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab 
%   License: BSD 3 clause

    if nargin < 14
        min_dist = 0.3;
        if nargin < 13
            verbose = false;
            if nargin < 12
                negative_sample_rate = 5;
                if nargin < 11
                    initial_alpha = 1;
                    if nargin < 10
                        gamma = 1;
                    end
                end
            end
        end
    end
    dim = size(head_embedding, 2);
    ONES=ones(1, dim);
    BG2S=2*gamma* b*ONES;
    FOURS=4*ONES;
    ABNEG2=-2.0*a*b;
    BNEG1=b-1;
    
    same_embedding = isequal(head_embedding, tail_embedding);
    alpha = initial_alpha;
    epochs_per_negative_sample = epochs_per_sample / single(negative_sample_rate);
    epoch_of_next_negative_sample = epochs_per_negative_sample;
    epoch_of_next_sample = epochs_per_sample;
    N=size(epochs_per_sample, 1);
    weights = ones(N,1)./epochs_per_sample; %We probably should have passed in weights to this instead...
    %tc=tic;
    
    TEST_CROSS_ENTROPY=false;
    if TEST_CROSS_ENTROPY
            dists = sqrt(sum((head_embedding(head,:) - tail_embedding(tail,:)).^2, 2));
            %spread = 1; %Should pass in spread from UMAP object.

            result = cross_entropy(double(head_embedding), double(tail_embedding), double(head), double(tail), double(weights), a, b, same_embedding);

            disp(['The current cross entropy is ' num2str(result)]);
    end
    for n = 1:n_epochs
        l=epoch_of_next_sample <= n;
        idxs=int32(find(l));
        N=int32(sum(l));
        nNegSamples = floor((single(n) - epoch_of_next_negative_sample(l)) ./ epochs_per_negative_sample(l));
        mxNeg=max(nNegSamples);
        randy=int32(randi(n_vertices, N, mxNeg));
        next=epoch_of_next_sample(l) + epochs_per_sample(l);
        next_neg=epoch_of_next_negative_sample(l)+(nNegSamples .* epochs_per_negative_sample(l));        
        %CANNOT vectorize dist_squared OUTSIDE next for loop since head_embedding and tail_embedding
        %change WITHIN next for loop since head and tail contain repeating values
        %   dist_squareds=sqrt(sum((head_embedding(head(l),:)-tail_embedding(tail(l),:)).^2, 2)).^2;
            
        for lIdx = 1:N
            i=idxs(lIdx);
            
            j = head(i);
            k = tail(i);
            
            current = head_embedding(j,:);
            other = tail_embedding(k,:);
            
            dist_squared = norm(current - other).^2;
            %CANNOT vectorize dist_squared OUTSIDE this loop
            %  assert(dist_squareds(li)==dist_squared);
            grad_coeff=(dist_squared > 0)*(ABNEG2*(dist_squared).^BNEG1)./ (a*dist_squared.^b + 1);
            grad_coeff(isnan(grad_coeff)) = 0;
            
            grad= max(-4, min(4, grad_coeff .* (current - other)));
            current = current + grad * alpha;
            if same_embedding
                other = other - grad * alpha;
                head_embedding(k,:) = other;
                tail_embedding(k,:) = other;
            end
                        
            n_neg_samples=nNegSamples(lIdx);
            for p = 1:n_neg_samples
                k=randy(lIdx,p);
                if same_embedding && j == k
                    continue
                end
                
                other = tail_embedding(k,:);
                
                dist_squared = norm(current - other).^2;
              
                grad_coeff=(dist_squared > 0)*BG2S./(0.001 + dist_squared)./(a*dist_squared.^b + 1);
                grad_coeff(isnan(grad_coeff)) = 0;
               
                grad=(grad_coeff > 0).*max(-4, min(4, grad_coeff .* (current - other))) + ~(grad_coeff > 0).*FOURS;
                
                current = current + grad * alpha;
            end
            
            head_embedding(j,:) = current;
            if same_embedding
                tail_embedding(j,:) = current;
            end
        end
        epoch_of_next_sample(l)=next;
        epoch_of_next_negative_sample(l)=next_neg;
        alpha = initial_alpha * (1 - single(n)/single(n_epochs));
        
        progress_checkpoint = min(floor(n_epochs / 10), 50);
        if verbose && mod(n, progress_checkpoint) == 0
            fprintf('\t%d/%d epochs done\n', int32(n), int32(n_epochs));
        end
    end
    
    embedding = head_embedding;
end    
    
