function ydata = tsne_bo(P, labels, no_dims)
%TSNE_P Performs symmetric t-SNE on affinity matrix P
if ~exist('labels', 'var')
    labels = [];
end
if ~exist('no_dims', 'var') || isempty(no_dims)
    no_dims = 2;
end

% First check whether we already have an initial solution
if numel(no_dims) > 1
    initial_solution = true;
    ydata = no_dims;
    no_dims = size(ydata, 2);
else
    initial_solution = false;
end
P = P / sum(sum(P));
% Initialize some variables
n = size(P, 1); % number of instances
momentum = .08; % initial momentum
final_momentum = .1; % value to which momentum is changed
mom_switch_iter = 250; % iteration at which momentum is changed
stop_lying_iter = 100; % iteration at which lying about P-values is stopped
max_iter = 1000; % maximum number of iterations
epsilon = 500; % initial learning rate
min_gain = .01; % minimum gain for delta-bar-delta

% Make sure P-vals are set properly
P(1:n+1:end) = 0; % set diagonal to zero
P = 0.5 * (P + P'); % symmetrize P-values
P = max(P./sum(P(:)), realmin); % make sure P-values sum to one
const = sum(P(:).*log(P(:))); % constant in KL divergence
if ~initial_solution
    P = P * 4; % lie about the P-vals to find better local minima
end

% Initialize the solution
if ~initial_solution
    ydata = .0001 * randn(n, no_dims);
end
y_incs = zeros(size(ydata));
gains = ones(size(ydata));

% Run the iterations
for iter = 1:max_iter

    % Compute joint probability that point i and j are neighbors
    sum_ydata = sum(ydata.^2, 2);
    num = 1 ./ (1 + bsxfun(@plus, sum_ydata, bsxfun(@plus, sum_ydata', -2*(ydata * ydata')))); % Student-t distribution
    num(1:n+1:end) = 0; % set diagonal to zero
    Q = max(num./sum(num(:)), realmin); % normalize to get probabilities

    % Compute the gradients (faster implementation)
    L = (P - Q) .* num;
    y_grads = 4 * (diag(sum(L, 1)) - L) * ydata;

    % Update the solution
    gains = (gains + .2) .* (sign(y_grads) ~= sign(y_incs)) ... % note that the y_grads are actually -y_grads
    +(gains * .8) .* (sign(y_grads) == sign(y_incs));
    gains(gains < min_gain) = min_gain;
    y_incs = momentum * y_incs - epsilon * (gains .* y_grads);
    ydata = ydata + y_incs;
    ydata = bsxfun(@minus, ydata, mean(ydata, 1));
    ydata(ydata < -100) = -100;
    ydata(ydata > 100) = 100;

    % Update the momentum if necessary
    if iter == mom_switch_iter
        momentum = final_momentum;
    end
    if iter == stop_lying_iter && ~initial_solution
        P = P ./ 4;
    end

    % Print out progress
    if ~rem(iter, 10)
        cost = const - sum(P(:).*log(Q(:)));
        %disp(['Iteration ' num2str(iter) ': error is ' num2str(cost)]);
    end

    % Display scatter plot (maximally first three dimensions)
    if ~isempty(labels)
        if no_dims == 1
            scatter(ydata, ydata, 9, labels, 'filled');
        elseif no_dims == 2
            scatter(ydata(:, 1), ydata(:, 2), 9, labels, 'filled');
        else
            scatter3(ydata(:, 1), ydata(:, 2), ydata(:, 3), 40, labels, 'filled');
        end
        axis equal tight
        %             axis off
        drawnow
    end
end
end