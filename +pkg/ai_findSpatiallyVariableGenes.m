function results = ai_findSpatiallyVariableGenes(expr, coords, nPerm)
% Identify spatially variable genes using Moran's I
% expr: genes x spots expression matrix
% coords: spots x 2 matrix of spatial coordinates (e.g., [x y])
% nPerm: number of permutations for significance testing

% https://chatgpt.com/share/6872a6bf-3c48-8005-8838-46e2f8a441e7

if nargin < 3
    nPerm = 1000; % default number of permutations
end

nGenes = size(expr,1);
nSpots = size(expr,2);

% Compute spatial weight matrix (inverse distance, zero diagonal)
D = squareform(pdist(coords));
W = 1 ./ (D + eps); % avoid division by zero
W(1:nSpots+1:end) = 0; % set diagonal to zero

results = table('Size', [nGenes, 3], ...
    'VariableTypes', {'double', 'double', 'double'}, ...
    'VariableNames', {'MoranI', 'pValue', 'qValue'});

for g = 1:nGenes
    x = expr(g,:)';
    x = x - mean(x);
    
    % Moran's I numerator and denominator
    num = x' * W * x;
    denom = x' * x;
    
    I = (nSpots / sum(W(:))) * (num / denom);
    results.MoranI(g) = I;
    
    % Permutation test
    permI = zeros(nPerm,1);
    for p = 1:nPerm
        xp = x(randperm(nSpots));
        permI(p) = (nSpots / sum(W(:))) * ((xp' * W * xp) / (xp' * xp));
    end
    
    results.pValue(g) = mean(permI >= I);
end

% Multiple testing correction
results.qValue = mafdr(results.pValue, 'BHFDR', true);
end


% expr: genes x spots matrix
% coords: spots x 2 spatial coordinates

% Example:

%{
expr = rand(1000, 250); % 1000 genes, 250 spots
coords = rand(250,2) * 100;

results = pkg.ai_findSpatiallyVariableGenes(expr, coords, 1000);
%}
% Spatially variable genes: e.g., results.qValue < 0.05


function results = findSpatiallyVariableGenesGP(expr, coords)
% Identify spatially variable genes using a simple GP test
% expr: genes x spots matrix
% coords: spots x 2 spatial coordinates

nGenes = size(expr,1);
nSpots = size(expr,2);

% Compute squared distance matrix
D = squareform(pdist(coords)).^2;

% Gaussian kernel bandwidth
sigma = median(D(:)); % heuristic: median distance squared
K = exp(-D / (2*sigma));

results = table('Size', [nGenes, 3], ...
    'VariableTypes', {'double', 'double', 'double'}, ...
    'VariableNames', {'LRTstat', 'pValue', 'qValue'});

for g = 1:nGenes
    y = expr(g,:)';
    y = y - mean(y); % center
    
    % Null model: noise only
    sigma2_null = var(y);
    LL_null = -0.5 * nSpots * log(2*pi*sigma2_null) - 0.5 * (y'*y) / sigma2_null;
    
    % GP model: spatial kernel + noise
    % variance parameters: sigma^2_gp and sigma^2_noise
    % Here use fixed sigma^2_gp = var(y)/2, sigma^2_noise = var(y)/2
    sigma2_gp = var(y)/2;
    sigma2_noise = var(y)/2;
    
    C = sigma2_gp * K + sigma2_noise * eye(nSpots);
    
    % Compute log likelihood
    L = chol(C + 1e-6*eye(nSpots),'lower'); % for numerical stability
    alpha = L'\(L\y);
    LL_gp = -0.5*(y'*alpha) - sum(log(diag(L))) - 0.5*nSpots*log(2*pi);
    
    % Likelihood ratio statistic
    LRT = 2*(LL_gp - LL_null);
    results.LRTstat(g) = LRT;
    
    % p-value: chi-squared with 1 df
    results.pValue(g) = 1 - chi2cdf(LRT,1);
end

% Adjust p-values
results.qValue = mafdr(results.pValue, 'BHFDR', true);
end



function results = runSpatialDE(expr, genes, coords)
% expr: genes x spots matrix
% genes: cell array of gene names
% coords: spots x 2 matrix

nGenes = size(expr,1);
nSpots = size(expr,2);

% Create sample names: spot_1, spot_2, ...
samples = strcat("spot_", string(1:nSpots))';

% Convert to long format
sample_col = repmat(samples, nGenes,1);
x_col = repmat(coords(:,1), nGenes,1);
y_col = repmat(coords(:,2), nGenes,1);
gene_col = repelem(string(genes), nSpots,1);
expr_col = reshape(expr', [],1);

% Create pandas DataFrame
pandas = py.importlib.import_module('pandas');
data = table(sample_col, x_col, y_col, gene_col, expr_col, ...
    'VariableNames', {'sample','x','y','gene','expr'});
df = pandas.DataFrame(data);

% Import SpatialDE
SpatialDE = py.importlib.import_module('SpatialDE');

% Step 1: preprocess (computes residuals)
resid = SpatialDE.preprocess(df);

% Step 2: run SpatialDE
results_df = SpatialDE.run(resid, df);

% Convert results back to MATLAB table
results = pandas2table(results_df);
end

function T = pandas2table(pdf)
% Convert pandas DataFrame to MATLAB table
cols = cell(pdf.columns.tolist());
data = cell(1,numel(cols));
for k = 1:numel(cols)
    colname = string(cols{k});
    data{k} = double(py.array.array('d', py.numpy.nditer(pdf{colname}.values)));
end
T = table(data{:}, 'VariableNames', matlab.lang.makeValidName(cols));
end
