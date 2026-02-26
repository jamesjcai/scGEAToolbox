function [Y] = sc_ifft(X, num_CCs_to_modify)
%scGFT — Synthetic Cell Generation
% https://github.com/Sanofi-Public/PMCB-scGFT/blob/master/R/Core.R

% https://pubmed.ncbi.nlm.nih.gov/39843603/
% https://pubmed.ncbi.nlm.nih.gov/31919373/
% https://pubmed.ncbi.nlm.nih.gov/37699885/
% https://pubmed.ncbi.nlm.nih.gov/35510186/

if nargin < 2, num_CCs_to_modify = 10; end

    [G, numc] = size(X);
    % G = number of genes (also components)
    % numc = number of cells

    if issparse(X), X = full(X); end
    X = pkg.norm_libsize(X, 1e4);
    Xn = log1p(X);            % X needs to be transposed next step!

    ft_mtx = fft(Xn.', [], 2);   % DFT along genes dimension, per cell, dim=2 is by rows
    % ft_mtx stores FFT freq. components
    assert(size(ft_mtx, 1) == numc)

    a = abs(ft_mtx);            % amplitudes
    a = a./vecnorm(a, 2, 1);    % Scaled amplitudes (Equation 20)
    sigma = var(a, 0, 1);

    synthesized_data = zeros(numc, G);

    valid_k = 1:floor((G-1)/2);
    num_to_mod = min(num_CCs_to_modify, length(valid_k));
    parent_idx = randi([1 numc], numc, 1);
    
    for cn = 1:numc
        cell_idx = parent_idx(cn);
        X_modified = ft_mtx(cn, :);
        selected_k = valid_k(randperm(length(valid_k), num_to_mod));
        for k = selected_k
            if size(synthesized_data, 1) > 1
                A_k = 2 + sqrt(sigma(k+1)) * randn();
            else
                A_k = 2 + randn();
            end
            
            % Apply modification to k-th component and its conjugate pair
            % X'[k] = A_k * X[k] (Equation 23)
            % Note: MATLAB uses 1-based indexing, so k+1 for component k
            X_modified(k + 1) = A_k * X_modified(k + 1);
            
            % Conjugate pair at N-k (which is G-k in 1-based indexing)
            % X'[N-k] = A_k * X[N-k]
            conj_idx = G - k + 1;
            X_modified(conj_idx) = A_k * X_modified(conj_idx);
        end        
        % Y(cn, :) = real(ifft(ft_mtx(cn, :)));
        synthesized_data(cn, :) = real(ifft(X_modified)); 
    end
    synthesized_data = max(0, synthesized_data);
    % Relative Percent Difference (RPD)
    % devs = mean(abs(Y - X));

    Y = convert_to_raw_counts(synthesized_data, X.', Xn.', parent_idx)';

    % if issparse(X)
    %     Y = sparse(cast(full(Y), 'like', X));
    % else
    %     Y = cast(Y, 'like', X);
    % end
end




function raw_synth = convert_to_raw_counts(synth_norm, raw_orig, norm_orig, parent_idx)
% Convert synthesized normalized data back to raw counts
% Uses relative change rate between synthesized and original normalized data
    
    num_synth = size(synth_norm, 1);
    num_genes = size(synth_norm, 2);
    raw_synth = zeros(num_synth, num_genes);
    
    for i = 1:num_synth
        orig_idx = parent_idx(i);
        
        % Get original normalized and raw values
        orig_norm = norm_orig(orig_idx, :);
        orig_raw = raw_orig(orig_idx, :);
        
        % Calculate relative change
        % Avoid division by zero
        orig_norm_safe = orig_norm;
        orig_norm_safe(orig_norm_safe == 0) = 1e-10;
        
        relative_change = synth_norm(i, :) ./ orig_norm_safe;
        
        % Apply to raw counts
        raw_synth(i, :) = round(orig_raw .* relative_change);
        
        % Handle cases where original was zero
        zero_mask = (orig_norm == 0);
        raw_synth(i, zero_mask) = 0;
    end    
    % Ensure non-negative integer counts
    raw_synth = max(0, round(raw_synth));
end


%% Wrapper function for single-cell synthesis
function [x_synth] = synthesize_single_cell(x_original, num_CCs)
% Synthesize a single cell from one original cell profile
% Simplified version for M=1 case
    
    G = length(x_original);
    
    % DFT
    X = fft(x_original);
    X_modified = X;
    
    % Identify valid components for modification
    if mod(G, 2) == 0
        valid_k = 1:(G/2 - 1);
    else
        valid_k = 1:floor((G-1)/2);
    end
    
    % Select and modify components
    num_to_mod = min(num_CCs, length(valid_k));
    selected_k = valid_k(randperm(length(valid_k), num_to_mod));
    
    for k = selected_k
        % For M=1: A_k ~ N(2, 1)
        A_k = 2 + randn();        
        X_modified(k + 1) = A_k * X_modified(k + 1);
        conj_idx = G - k + 1;
        X_modified(conj_idx) = A_k * X_modified(conj_idx);
    end
    
    % IFFT and ReLU
    x_synth = real(ifft(X_modified));
    x_synth = max(0, x_synth);
end


%{
https://chatgpt.com/s/t_699b270e1f2081919a22b558d6bb6e47

scGFT Synthetic Cell Generation

Input dataset: Lung_dev_scRNA (12,384 cells)

Generate:
  (•) 50% additional cells
  ( ) Custom number: [      ]

Mode:
  (•) Per cluster
  ( ) Global

Rare cell amplification:
  [x] Oversample clusters with < 2% cells

Output:
  (•) Append to dataset
  ( ) Create new dataset
%}