function [activities, stats] = ai_ULM(expressionMatrix, regulons, varargin)
% ULM - Univariate Linear Model for regulator activity inference
%
% Implementation of the ULM method from decoupleR for inferring transcription
% factor or pathway activities from gene expression data using univariate
% linear models.
%
% INPUTS:
%   expressionMatrix - genes x samples matrix of expression values
%   regulons - struct array with fields:
%              'name' - regulator name
%              'targets' - vector of target gene indices
%              'weights' - optional weights for targets (default: all 1)
%   Optional parameters:
%     'geneNames' - cell array of gene names
%     'sampleNames' - cell array of sample names
%     'minTargets' - minimum number of targets required (default: 5)
%     'alpha' - significance level for p-values (default: 0.05)
%     'verbose' - logical, display progress (default: true)
%     'standardize' - logical, standardize expression (default: true)
%
% OUTPUTS:
%   activities - struct with fields:
%                'scores' - regulators x samples matrix of activity scores
%                'pvalues' - regulators x samples matrix of p-values
%                'regulatorNames' - cell array of regulator names
%                'sampleNames' - cell array of sample names
%   stats - struct with additional statistics

% Parse inputs
p = inputParser;
addRequired(p, 'expressionMatrix', @isnumeric);
addRequired(p, 'regulons', @isstruct);
addParameter(p, 'geneNames', {}, @iscell);
addParameter(p, 'sampleNames', {}, @iscell);
addParameter(p, 'minTargets', 5, @isnumeric);
addParameter(p, 'alpha', 0.05, @isnumeric);
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'standardize', true, @islogical);
parse(p, expressionMatrix, regulons, varargin{:});

[nGenes, nSamples] = size(expressionMatrix);
nRegulators = length(regulons);

if p.Results.verbose
    fprintf('ULM Analysis:\n');
    fprintf('  Genes: %d\n', nGenes);
    fprintf('  Samples: %d\n', nSamples);
    fprintf('  Regulators: %d\n', nRegulators);
end

% Standardize expression matrix if requested
if p.Results.standardize
    expressionMatrix = zscore(expressionMatrix, 0, 2); % Standardize across samples
end

% Initialize output matrices
activityScores = zeros(nRegulators, nSamples);
pValues = zeros(nRegulators, nSamples);
validRegulators = true(nRegulators, 1);
nTargetsUsed = zeros(nRegulators, 1);

% Create regulator names
if isfield(regulons, 'name')
    regulatorNames = {regulons.name};
else
    regulatorNames = arrayfun(@(x) sprintf('Regulator_%d', x), 1:nRegulators, 'UniformOutput', false);
end

% Create sample names if not provided
if isempty(p.Results.sampleNames)
    sampleNames = arrayfun(@(x) sprintf('Sample_%d', x), 1:nSamples, 'UniformOutput', false);
else
    sampleNames = p.Results.sampleNames;
end

% Process each regulator
if p.Results.verbose
    fprintf('Processing regulators...\n');
end

for regIdx = 1:nRegulators
    if p.Results.verbose && mod(regIdx, 50) == 0
        fprintf('  Processing regulator %d/%d\n', regIdx, nRegulators);
    end
    
    targets = regulons(regIdx).targets;
    
    % Remove invalid targets
    targets = targets(targets >= 1 & targets <= nGenes);
    
    % Check minimum targets requirement
    if length(targets) < p.Results.minTargets
        if p.Results.verbose
            fprintf('    Skipping regulator %s: only %d targets (min: %d)\n', ...
                regulatorNames{regIdx}, length(targets), p.Results.minTargets);
        end
        validRegulators(regIdx) = false;
        continue;
    end
    
    nTargetsUsed(regIdx) = length(targets);
    
    % Get weights if available
    if isfield(regulons, 'weights') && ~isempty(regulons(regIdx).weights)
        weights = regulons(regIdx).weights;
        weights = weights(1:length(targets)); % Ensure same length as targets
    else
        weights = ones(length(targets), 1);
    end
    
    % Get target expression
    targetExpression = expressionMatrix(targets, :);
    
    % Apply ULM for this regulator
    [scores, pvals] = fitULM(targetExpression, weights);
    
    activityScores(regIdx, :) = scores;
    pValues(regIdx, :) = pvals;
end

% Filter out invalid regulators
activityScores = activityScores(validRegulators, :);
pValues = pValues(validRegulators, :);
regulatorNames = regulatorNames(validRegulators);
nTargetsUsed = nTargetsUsed(validRegulators);

% Prepare output
activities = struct();
activities.scores = activityScores;
activities.pvalues = pValues;
activities.regulatorNames = regulatorNames;
activities.sampleNames = sampleNames;

% Additional statistics
stats = struct();
stats.nValidRegulators = sum(validRegulators);
stats.nTargetsUsed = nTargetsUsed;
stats.minTargets = p.Results.minTargets;
stats.alpha = p.Results.alpha;

% Calculate summary statistics
stats.meanActivity = mean(activityScores, 2);
stats.stdActivity = std(activityScores, 0, 2);
stats.significantSamples = sum(pValues <= p.Results.alpha, 2);

if p.Results.verbose
    fprintf('ULM analysis completed:\n');
    fprintf('  Valid regulators: %d/%d\n', stats.nValidRegulators, nRegulators);
    fprintf('  Mean targets per regulator: %.1f\n', mean(nTargetsUsed));
end

end

function [activityScores, pValues] = fitULM(targetExpression, weights)
% Fit Univariate Linear Model for a single regulator
%
% INPUTS:
%   targetExpression - targets x samples matrix
%   weights - vector of target weights
%
% OUTPUTS:
%   activityScores - vector of activity scores across samples
%   pValues - vector of p-values for each sample

[nTargets, nSamples] = size(targetExpression);

% Normalize weights
weights = weights / sum(weights);
weights = weights(:); % Ensure column vector

activityScores = zeros(1, nSamples);
pValues = zeros(1, nSamples);

% For each sample, estimate activity using weighted approach
for sampleIdx = 1:nSamples
    geneExpr = targetExpression(:, sampleIdx);
    
    % Calculate weighted activity score
    activityScore = sum(weights .* geneExpr);
    
    % Calculate statistics for significance testing
    % Use one-sample t-test against zero
    if nTargets > 1
        % Weighted standard error
        weightedMean = activityScore;
        weightedVar = sum(weights.^2 .* (geneExpr - weightedMean).^2);
        se = sqrt(weightedVar);
        
        if se > 0
            tStat = weightedMean / se;
            df = nTargets - 1;
            pValue = 2 * (1 - tcdf(abs(tStat), df)); % Two-tailed test
        else
            pValue = 1;
        end
    else
        pValue = 1;
    end
    
    activityScores(sampleIdx) = activityScore;
    pValues(sampleIdx) = pValue;
end

end

function regulons = createRegulon(regulatorName, targetGenes, weights)
% Helper function to create a single regulon structure
%
% INPUTS:
%   regulatorName - string, name of the regulator
%   targetGenes - vector of target gene indices
%   weights - optional vector of weights (default: all 1s)

if nargin < 3
    weights = ones(length(targetGenes), 1);
end

regulons = struct();
regulons.name = regulatorName;
regulons.targets = targetGenes(:);
regulons.weights = weights(:);

end

function regulons = loadRegulonsFromTable(regulonTable, allGeneNames)
% Convert regulon table to regulon structure
%
% INPUTS:
%   regulonTable - table with columns 'regulator', 'target', 'weight' (optional)
%   allGeneNames - cell array of all gene names
%
% OUTPUT:
%   regulons - struct array of regulons

if ~istable(regulonTable)
    error('regulonTable must be a MATLAB table');
end

% Get unique regulators
uniqueRegulators = unique(regulonTable.regulator);
nRegulators = length(uniqueRegulators);

regulons = struct('name', {}, 'targets', {}, 'weights', {});

for i = 1:nRegulators
    regName = uniqueRegulators{i};
    
    % Get targets for this regulator
    regRows = strcmp(regulonTable.regulator, regName);
    targetNames = regulonTable.target(regRows);
    
    % Convert gene names to indices
    [~, targetIndices] = ismember(targetNames, allGeneNames);
    validTargets = targetIndices > 0;
    
    if sum(validTargets) == 0
        continue;
    end
    
    targetIndices = targetIndices(validTargets);
    
    % Get weights if available
    if ismember('weight', regulonTable.Properties.VariableNames)
        targetWeights = regulonTable.weight(regRows);
        targetWeights = targetWeights(validTargets);
    else
        targetWeights = ones(sum(validTargets), 1);
    end
    
    % Create regulon
    regulons(end+1) = createRegulon(regName, targetIndices, targetWeights);
end

end

