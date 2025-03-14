function [Y] = sc_diffumap(X)
    %X = sc_transform(X);
    X = X.';
    % n-p
    [Y] = standardDiffusionMap(X, 0.05, 100, 2);
end

% Implementation of both Standard and Fokker-Planck Diffusion Maps

function [phi, lambda] = standardDiffusionMap(X, epsilon, t, n_components)
    % Standard Diffusion Map implementation
    % X: data matrix (n_samples x n_features)
    % epsilon: kernel bandwidth
    % t: diffusion time
    % n_components: number of components to return
    
    % Compute pairwise distances
    D = pdist2(X, X, 'cosine');
    
    % Construct kernel matrix
    K = exp(-D.^2 / (2*epsilon));
    
    % Compute row sums for normalization
    d = sum(K, 2);
    
    % Normalize to create transition matrix
    P = K ./ (d * d');
    P = P ./ sum(P, 2);
    
    % Compute eigenvalues and eigenvectors
    [V, D] = eigs(P, n_components+1);
    
    % Sort eigenvalues and eigenvectors
    [lambda, idx] = sort(diag(D), 'descend');
    V = V(:, idx);
    
    % Remove first trivial eigenvector and eigenvalue
    lambda = lambda(2:end);
    phi = V(:, 2:end);
    
    % Scale eigenvectors by eigenvalues
    phi = phi * diag(lambda.^t);
end

%{

function [phi, lambda] = fokkerPlanckDiffusionMap(X, epsilon, dt, n_components)
    % Fokker-Planck Diffusion Map implementation
    % X: data matrix (n_samples x n_features)
    % epsilon: kernel bandwidth
    % dt: time step
    % n_components: number of components to return
    
    % Compute pairwise distances
    D = pdist2(X, X, 'euclidean');
    
    % Estimate drift term (using finite differences)
    n_samples = size(X, 1);
    drift = zeros(size(X));
    for i = 1:n_samples
        % Find k nearest neighbors
        [~, idx] = sort(D(i,:));
        k = min(10, n_samples-1); % Use 10 nearest neighbors
        neighbors = idx(2:k+1);
        
        % Estimate drift using local averaging
        drift(i,:) = mean(X(neighbors,:) - X(i,:), 1) / dt;
    end
    
    % Construct asymmetric kernel incorporating drift
    K = zeros(n_samples);
    for i = 1:n_samples
        for j = 1:n_samples
            % Add drift contribution to standard kernel
            drift_term = dot(drift(i,:), X(j,:) - X(i,:));
            K(i,j) = exp(-D(i,j)^2/(2*epsilon) + drift_term*dt);
        end
    end
    
    % Normalize to create forward operator
    d = sum(K, 2);
    P = K ./ (d * d');
    P = P ./ sum(P, 2);
    
    % Compute eigenvalues and eigenvectors
    [V, D] = eigs(P, n_components+1);
    
    % Sort eigenvalues and eigenvectors
    [lambda, idx] = sort(diag(D), 'descend');
    V = V(:, idx);
    
    % Remove first trivial eigenvector and eigenvalue
    lambda = lambda(2:end);
    phi = V(:, 2:end);
    
    % Scale eigenvectors by eigenvalues
    phi = phi * diag(lambda.^(1/dt));
end

%}
% Example usage:
% X = randn(100, 3);  % Generate random 3D data
% epsilon = 1;        % Kernel bandwidth
% t = 1;             % Diffusion time
% n_components = 2;   % Number of components to keep
%
% % Compute standard diffusion map
% [phi_std, lambda_std] = standardDiffusionMap(X, epsilon, t, n_components);
%
% % Compute Fokker-Planck diffusion map
% [phi_fp, lambda_fp] = fokkerPlanckDiffusionMap(X, epsilon, t, n_components);


%{
No, a Fokker-Planck diffusion map is different from a standard diffusion map. While both methods are related to analyzing dynamics and probability distributions, they have distinct characteristics:
Standard diffusion maps (Coifman & Lafon) focus on:

Constructing similarity matrices based on local geometric structure
Using eigenvectors of the normalized graph Laplacian
Primarily used for dimensionality reduction and data visualization

Fokker-Planck diffusion maps involve:

Incorporating drift terms from the Fokker-Planck equation
Handling non-equilibrium dynamics
Can capture asymmetric transitions and directed flows

The key difference is that standard diffusion maps assume symmetric diffusion processes, while Fokker-Planck diffusion maps can handle asymmetric transitions and non-equilibrium dynamics. Would you like me to elaborate on any particular aspect of these differences?

```matlab
% Implementation of both Standard and Fokker-Planck Diffusion Maps

function [phi, lambda] = standardDiffusionMap(X, epsilon, t, n_components)
    % Standard Diffusion Map implementation
    % X: data matrix (n_samples x n_features)
    % epsilon: kernel bandwidth
    % t: diffusion time
    % n_components: number of components to return
    
    % Compute pairwise distances
    D = pdist2(X, X, 'euclidean');
    
    % Construct kernel matrix
    K = exp(-D.^2 / (2*epsilon));
    
    % Compute row sums for normalization
    d = sum(K, 2);
    
    % Normalize to create transition matrix
    P = K ./ (d * d');
    P = P ./ sum(P, 2);
    
    % Compute eigenvalues and eigenvectors
    [V, D] = eigs(P, n_components+1);
    
    % Sort eigenvalues and eigenvectors
    [lambda, idx] = sort(diag(D), 'descend');
    V = V(:, idx);
    
    % Remove first trivial eigenvector and eigenvalue
    lambda = lambda(2:end);
    phi = V(:, 2:end);
    
    % Scale eigenvectors by eigenvalues
    phi = phi * diag(lambda.^t);
end

function [phi, lambda] = fokkerPlanckDiffusionMap(X, epsilon, dt, n_components)
    % Fokker-Planck Diffusion Map implementation
    % X: data matrix (n_samples x n_features)
    % epsilon: kernel bandwidth
    % dt: time step
    % n_components: number of components to return
    
    % Compute pairwise distances
    D = pdist2(X, X, 'euclidean');
    
    % Estimate drift term (using finite differences)
    n_samples = size(X, 1);
    drift = zeros(size(X));
    for i = 1:n_samples
        % Find k nearest neighbors
        [~, idx] = sort(D(i,:));
        k = min(10, n_samples-1); % Use 10 nearest neighbors
        neighbors = idx(2:k+1);
        
        % Estimate drift using local averaging
        drift(i,:) = mean(X(neighbors,:) - X(i,:), 1) / dt;
    end
    
    % Construct asymmetric kernel incorporating drift
    K = zeros(n_samples);
    for i = 1:n_samples
        for j = 1:n_samples
            % Add drift contribution to standard kernel
            drift_term = dot(drift(i,:), X(j,:) - X(i,:));
            K(i,j) = exp(-D(i,j)^2/(2*epsilon) + drift_term*dt);
        end
    end
    
    % Normalize to create forward operator
    d = sum(K, 2);
    P = K ./ (d * d');
    P = P ./ sum(P, 2);
    
    % Compute eigenvalues and eigenvectors
    [V, D] = eigs(P, n_components+1);
    
    % Sort eigenvalues and eigenvectors
    [lambda, idx] = sort(diag(D), 'descend');
    V = V(:, idx);
    
    % Remove first trivial eigenvector and eigenvalue
    lambda = lambda(2:end);
    phi = V(:, 2:end);
    
    % Scale eigenvectors by eigenvalues
    phi = phi * diag(lambda.^(1/dt));
end

% Example usage:
% X = randn(100, 3);  % Generate random 3D data
% epsilon = 1;        % Kernel bandwidth
% t = 1;             % Diffusion time
% n_components = 2;   % Number of components to keep
%
% % Compute standard diffusion map
% [phi_std, lambda_std] = standardDiffusionMap(X, epsilon, t, n_components);
%
% % Compute Fokker-Planck diffusion map
% [phi_fp, lambda_fp] = fokkerPlanckDiffusionMap(X, epsilon, t, n_components);

```

I've implemented both standard diffusion maps and Fokker-Planck diffusion maps in MATLAB. Here are the key differences in the implementations:

1. Standard Diffusion Map:
- Uses symmetric kernel based on Euclidean distances
- Implements the classic Coifman & Lafon algorithm
- Simpler normalization procedure

2. Fokker-Planck Diffusion Map:
- Incorporates drift terms estimated from local neighborhoods
- Uses an asymmetric kernel that accounts for directional flow
- Includes time step parameter for evolution

The code includes example usage at the bottom. Would you like me to explain any particular part of the implementation in more detail?

Some potential modifications you might want to consider:
- Changing the kernel function (e.g., using different distance metrics)
- Modifying the drift estimation method
- Adding adaptive bandwidth selection
- Including different normalization schemes

%}