% https://chat.tamu.ai/s/3f32af28-cf33-4aeb-8a6d-c50f78e8dc18

% ---- Small synthetic example (G genes, C cells) ----
rng(1);

G = 50;
C = 200;

% Create 3 "true" cellstates with different gene fractions
Ktrue = 3;
phiTrue = rand(G, Ktrue);
phiTrue = phiTrue ./ sum(phiTrue,1);

% Assign cells to true states
zTrue = randi(Ktrue, 1, C);

% Sample total UMIs per cell
Nc = poissrnd(1500, 1, C) + 200;   % around 1700 UMIs/cell

% Generate counts from multinomial
N = sparse(G, C);
for c = 1:C
    p = phiTrue(:, zTrue(c));
    N(:,c) = mnrnd(Nc(c), p)';     % Statistics Toolbox
end
N = sparse(N);

% ---- Run CELLSTATES partition optimization (fixed Theta) ----
%Theta = 2^11;  % example scale; paper uses Theta=2^q and starts near avg NUMI [2]
Theta = 4096;

opts = struct();
opts.S = C;        % default in supplement [2]
opts.T = 1000;     % default in supplement [2]
opts.verbose = true;
opts.doPost = true;

out = cellstates_partition_opt(N, Theta, opts);

% Results
labels = out.labels;   % 1..K inferred clusters
K = out.K;
fprintf('Inferred K = %d clusters, logL = %.3f\n', K, out.logL);

% Optional: compare to truth (quick purity-like check)
tab = crosstab(zTrue(:), labels(:));
disp(tab);