function [n,result] = howmany(alpha,R,J)
    % howmany computes the actuan simulated annealing (SA)
    % INPUT:
    % R =========> Redundancy matrix (MI of feature-feature
    % J =========> Importance vector (MI of target-features)
    % alpha =====> Mixing optimal parameter
    % OUTPUT
    % results ===> QUBO object with solution information
    % n =========> Number of selected features

    % Compute Q matrix
    Q = qubo((1-alpha)*R - alpha*diag(J));
    % QUBO object with SA solutions 
    result = solve(Q);
    % Number of features found
    n = sum(result.BestX);
end