function [aucScores, rankings] = ai_AUCell(expressionMatrix, geneSets, varargin)
% AUCell - Calculate Area Under the Curve scores for gene sets in single cells
%
% INPUTS:
%   expressionMatrix - genes x cells matrix of expression values
%   geneSets - cell array where each element is a vector of gene indices
%              or a struct with fields 'genes' (gene indices) and 'name'
%   Optional parameters:
%     'aucMaxRank' - maximum rank to consider for AUC calculation (default: 5% of genes)
%     'geneNames' - cell array of gene names (for validation)
%     'verbose' - logical, display progress (default: true)
%
% OUTPUTS:
%   aucScores - gene sets x cells matrix of AUC scores
%   rankings - genes x cells matrix of gene rankings per cell

% Parse inputs
p = inputParser;
addRequired(p, 'expressionMatrix', @isnumeric);
addRequired(p, 'geneSets');
addParameter(p, 'aucMaxRank', [], @isnumeric);
addParameter(p, 'geneNames', {}, @iscell);
addParameter(p, 'verbose', true, @islogical);
parse(p, expressionMatrix, geneSets, varargin{:});

% Get dimensions
[nGenes, nCells] = size(expressionMatrix);

% Set default aucMaxRank (5% of genes, minimum 50)
if isempty(p.Results.aucMaxRank)
    aucMaxRank = max(50, round(0.05 * nGenes));
else
    aucMaxRank = p.Results.aucMaxRank;
end

if p.Results.verbose
    fprintf('AUCell analysis:\n');
    fprintf('  Genes: %d\n', nGenes);
    fprintf('  Cells: %d\n', nCells);
    fprintf('  AUC max rank: %d\n', aucMaxRank);
end

% Step 1: Create rankings for each cell
if p.Results.verbose
    fprintf('Creating gene rankings...\n');
end

rankings = zeros(size(expressionMatrix));
for i = 1:nCells
    [~, rankings(:,i)] = sort(expressionMatrix(:,i), 'descend');
end

% Convert gene sets to consistent format
if isstruct(geneSets)
    geneSetList = geneSets;
    nGeneSets = length(geneSetList);
else
    % Convert cell array to struct format
    nGeneSets = length(geneSets);
    geneSetList = struct('genes', geneSets, 'name', []);
    for i = 1:nGeneSets
        if isempty(geneSetList(i).name)
            geneSetList(i).name = sprintf('GeneSet_%d', i);
        end
    end
end

% Initialize output
aucScores = zeros(nGeneSets, nCells);

% Step 2: Calculate AUC for each gene set in each cell
if p.Results.verbose
    fprintf('Calculating AUC scores...\n');
end

for setIdx = 1:nGeneSets
    geneSet = geneSetList(setIdx).genes;
    
    % Validate gene set
    if any(geneSet < 1 | geneSet > nGenes)
        warning('Gene set %d contains invalid gene indices', setIdx);
        continue;
    end
    
    if p.Results.verbose && mod(setIdx, 10) == 0
        fprintf('  Processing gene set %d/%d\n', setIdx, nGeneSets);
    end
    
    % Calculate AUC for this gene set across all cells
    for cellIdx = 1:nCells
        aucScores(setIdx, cellIdx) = calculateAUC(rankings(:, cellIdx), geneSet, aucMaxRank, nGenes);
    end
end

if p.Results.verbose
    fprintf('AUCell analysis complete!\n');
end

end

function auc = calculateAUC(geneRanks, geneSet, maxRank, nGenes)
% Calculate AUC for a single gene set in a single cell
%
% INPUTS:
%   geneRanks - vector of gene rankings (1 = highest expressed)
%   geneSet - vector of gene indices in the set
%   maxRank - maximum rank to consider
%   nGenes - total number of genes

% Get ranks of genes in the set
setRanks = geneRanks(geneSet);

% Remove genes not in top maxRank
setRanks = setRanks(setRanks <= maxRank);

if isempty(setRanks)
    auc = 0;
    return;
end

% Sort the ranks
setRanks = sort(setRanks);
nGenesInSet = length(setRanks);

% Calculate AUC using trapezoidal rule
% x-axis: rank position (1 to maxRank)
% y-axis: cumulative recovery (0 to 1)

x = [0; setRanks; maxRank];
y = [0; (1:nGenesInSet)'/nGenesInSet; 1];

% Calculate area under the curve
auc = trapz(x, y) / maxRank;

% Normalize: subtract the area under random curve (0.5)
% and scale to [0,1] range
auc = max(0, (auc - 0.5) * 2);

end

% Helper function to create gene sets from gene names
function geneSets = createGeneSetsFromNames(geneSetNames, allGeneNames)
% Convert gene sets defined by gene names to gene indices
%
% INPUTS:
%   geneSetNames - cell array where each element is a cell array of gene names
%   allGeneNames - cell array of all gene names in the expression matrix
%
% OUTPUT:
%   geneSets - cell array of gene indices

nSets = length(geneSetNames);
geneSets = cell(nSets, 1);

for i = 1:nSets
    currentSet = geneSetNames{i};
    [~, geneIndices] = ismember(currentSet, allGeneNames);
    
    % Remove genes not found
    geneIndices = geneIndices(geneIndices > 0);
    
    if isempty(geneIndices)
        warning('Gene set %d: No genes found in expression matrix', i);
    end
    
    geneSets{i} = geneIndices;
end

end

