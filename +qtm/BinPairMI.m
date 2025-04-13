function MI = BinPairMI(x, y)
    % BinPairMI computes the mutual information within x and y
    % INPUT:
    % x ====> Feature x observations
    % y ====> Feature y observations
    % OUTPUT
    % MI ===> Mutual information value within x-y

    % Estimate joint probability distribution with MATLAB automatic binning
    [joint_counts, ~, ~] = histcounts2(x, y);
  
    % Try to avoid singularities with this and eps0
    nconts = sum(joint_counts(:));
    if nconts == 0
        MI = 0;
        return;
    end

    % Handle zeros (add eps for numerical stability)
    eps0 = eps(realmin('single'));

    joint_prob = joint_counts ./ nconts + eps0;
   
    % Estimate marginal probabilities
    x_marginal = sum(joint_prob, 2) + eps0;
    y_marginal = sum(joint_prob, 1) + eps0;
    
    % Calculate entropy terms
    entropy_xy = -sum(joint_prob(:) .* log2(joint_prob(:)) );
    entropy_x = -sum(x_marginal .* log2(x_marginal) );
    entropy_y = -sum(y_marginal .* log2(y_marginal) );
    
    % Mutual information
    MI = entropy_x + entropy_y - entropy_xy;
end