function [X] = norm_gficf(X)

% ref: https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2019.00734/full
    % UMI cell count matrix gene-by-cell
    % X = cpm(calcNormFactors(X)); % Normalize counts using EdgeR-like normalization
    [X] = pkg.norm_libsize(X, 1e6);

    % Compute GF transformation
    X = X ./ sum(X, 1); % Normalize by column sums
    
    % Get IDF weights
    w = getIdfW(X);
    
    % Apply ICF
    X = X.*w;

    % Normalize with L2 norm
    % gficf = l_norm(gficf, 'l2');
    X = X./vecnorm(X, 1);
end


function w = getIdfW(M, type)
    nt = sum(M ~= 0, 2);
    switch type
        case 'classic'
            w = log((size(M, 2) + 1) ./ (nt + 1));
        case 'prob'
            w = log((size(M, 2) - nt) ./ nt);
        case 'smooth'
            w = log(1 + size(M, 2) ./ nt);
    end
end

% function M = l_norm(m, norm)
%     if strcmp(norm, 'l1')
%         norm_vec = 1 ./ sum(m, 2);
%     elseif strcmp(norm, 'l2')
%         norm_vec = 1 ./ sqrt(sum(m.^2, 2));
%     end
%     norm_vec(isinf(norm_vec)) = 0; % Handle infinite values
% 
%     % Apply normalization
%     M = m * diag(norm_vec); % Assuming sparse input
% end