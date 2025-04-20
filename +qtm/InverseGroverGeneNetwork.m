% Inverse Grover's Algorithm for Gene Expression Analysis
% This implementation discovers significant gene expression patterns from data

classdef InverseGroverGeneNetwork
    properties
        numGenes;           % Number of genes in the network
        expressionMatrix;   % Gene expression matrix (genes x samples)
        geneNames;          % Names of genes
        threshold;          % Threshold for binarizing expression values
        binaryMatrix;       % Binarized expression matrix
        uniquePatterns;     % Unique expression patterns found in data
    end
    
    methods
        function obj = InverseGroverGeneNetwork(expressionMatrix, geneNames, threshold)
            % Constructor for the InverseGroverGeneNetwork class
            %
            % Args:
            %   expressionMatrix: Gene expression matrix (genes x samples)
            %   geneNames: Cell array of gene names
            %   threshold: Threshold for binarizing expression (default: 0.5)
            
            if nargin < 3
                threshold = 0.5;
            end
            
            obj.expressionMatrix = expressionMatrix;
            obj.numGenes = size(expressionMatrix, 1);
            
            % Ensure geneNames is a cell array of the correct size
            if ischar(geneNames)
                obj.geneNames = {geneNames};
            elseif iscell(geneNames)
                obj.geneNames = geneNames;
            else
                obj.geneNames = arrayfun(@(i) sprintf('Gene%d', i), 1:obj.numGenes, 'UniformOutput', false);
            end
            
            % Make sure geneNames has the correct length
            if length(obj.geneNames) ~= obj.numGenes
                obj.geneNames = arrayfun(@(i) sprintf('Gene%d', i), 1:obj.numGenes, 'UniformOutput', false);
            end
            
            obj.threshold = threshold;
            
            % Binarize the expression matrix
            obj.binaryMatrix = double(obj.expressionMatrix > threshold);
            
            % Find unique patterns
            obj.uniquePatterns = unique(obj.binaryMatrix', 'rows');
        end
        
        function validationResult = validationFunction(obj, stateStr)
            % Validation function that checks if a state matches any pattern in the data
            %
            % Args:
            %   stateStr: String representation of a state (e.g., '101')
            %
            % Returns:
            %   validationResult: True if state matches a pattern, false otherwise
            
            % Convert state string to numeric array
            stateArray = zeros(1, obj.numGenes);
            for i = 1:obj.numGenes
                stateArray(i) = str2double(stateStr(i));
            end
            
            % Check if state matches any unique pattern
            validationResult = any(all(obj.uniquePatterns == stateArray, 2));
        end
        
        function U = createOracle(obj)
            % Creates an oracle matrix that applies a phase flip to marked states
            %
            % Returns:
            %   U: Unitary matrix representing the oracle
            
            % Create matrix of size 2^n x 2^n
            n = obj.numGenes;
            N = 2^n;
            U = eye(N);
            
            % Apply phase flip to states that match patterns in the data
            for i = 0:(N-1)
                % Convert decimal i to binary string
                binaryStr = dec2bin(i, n);
                
                % Apply validation function
                if obj.validationFunction(binaryStr)
                    % Apply phase flip by multiplying by -1
                    U(i+1, i+1) = -1;
                end
            end
        end
        
        function U = createDiffusion(obj)
            % Creates the diffusion operator matrix
            %
            % Returns:
            %   U: Unitary matrix representing the diffusion operator
            
            n = obj.numGenes;
            N = 2^n;
            
            % Create matrix for state |s>
            s = ones(N, 1) / sqrt(N);
            
            % U = 2|s><s| - I
            U = 2 * (s * s') - eye(N);
        end
        
        function [discoveredPatterns, probabilities] = runInverseGrover(obj, numIterations)
            % Runs inverse Grover's algorithm to discover patterns
            %
            % Args:
            %   numIterations: Number of Grover iterations (if not specified, 
            %                  will use optimal estimate)
            %
            % Returns:
            %   discoveredPatterns: Binary strings representing discovered patterns
            %   probabilities: Corresponding probability for each pattern
            
            n = obj.numGenes;
            N = 2^n;
            
            % Estimate number of marked states (unique patterns)
            k = size(obj.uniquePatterns, 1);
            k = max(1, k); % Ensure k is at least 1 to avoid division by zero
            
            % Calculate optimal number of iterations if not provided
            if nargin < 2 || isempty(numIterations)
                numIterations = floor(pi/4 * sqrt(N/k));
                numIterations = max(1, numIterations); % Ensure at least 1 iteration
            end
            
            % Create oracle and diffusion operators
            oracleU = obj.createOracle();
            diffusionU = obj.createDiffusion();
            
            % Initialize state vector with equal superposition
            psi = ones(N, 1) / sqrt(N);
            
            numIterations=3;
            
            % Apply Grover iterations
            for i = 1:numIterations
                % Apply oracle
                psi = oracleU * psi;
                
                % Apply diffusion
                psi = diffusionU * psi;
                
                % Optional: Display progress
                fprintf('Completed iteration %d of %d\n', i, numIterations);
            end
            
            % Calculate probabilities of measuring each state
            probArray = abs(psi).^2;
            
            % Get the most probable states (above threshold)
            threshold = 1/(2*N); % Probability threshold
            highProbIndices = find(probArray > threshold);
            
            % Convert indices to binary patterns
            discoveredPatterns = cell(length(highProbIndices), 1);
            probabilities = zeros(length(highProbIndices), 1);
            
            for i = 1:length(highProbIndices)
                idx = highProbIndices(i) - 1; % Convert to 0-based index
                discoveredPatterns{i} = dec2bin(idx, n);
                probabilities(i) = probArray(highProbIndices(i));
            end
            
            % Sort by probability (descending)
          %  [probabilities, sortIdx] = sort(probabilities, 'descend');
          %  discoveredPatterns = discoveredPatterns(sortIdx);
        end
        
        function results = analyzePatterns(obj, discoveredPatterns, probabilities)
            % Analyzes discovered patterns and maps them to gene information
            %
            % Args:
            %   discoveredPatterns: Cell array of discovered pattern strings
            %   probabilities: Corresponding probabilities
            %
            % Returns:
            %   results: Structure with analysis results
            
            numPatterns = length(discoveredPatterns);
            patternAnalysis = cell(numPatterns, 1);
            
            for i = 1:numPatterns
                pattern = discoveredPatterns{i};
                expressed = {};
                suppressed = {};
                
                % Identify expressed and suppressed genes
                for j = 1:obj.numGenes
                    if pattern(j) == '1'
                        expressed{end+1} = obj.geneNames{j};
                    else
                        suppressed{end+1} = obj.geneNames{j};
                    end
                end
                
                % Create analysis struct for this pattern
                analysisStruct = struct;
                analysisStruct.pattern = pattern;
                analysisStruct.probability = probabilities(i);
                analysisStruct.expressedGenes = expressed;
                analysisStruct.suppressedGenes = suppressed;
                
                % Store analysis for this pattern
                patternAnalysis{i} = analysisStruct;
            end
            
            % Create final results structure
            results = struct;
            results.numGenesAnalyzed = obj.numGenes;
            results.numSamples = size(obj.expressionMatrix, 2);
            results.numPatternsDiscovered = numPatterns;
            results.patternAnalysis = patternAnalysis;
        end
        
        function visualizeResults(obj, results)
            % Visualizes the results of pattern discovery
            %
            % Args:
            %   results: Results structure from analyzePatterns
            
            numPatterns = results.numPatternsDiscovered;
            
            if numPatterns == 0
                fprintf('No significant patterns discovered.\n');
                return;
            end
            
            patterns = cell(numPatterns, 1);
            probs = zeros(numPatterns, 1);
            
            % Extract patterns and probabilities for plotting
            for i = 1:numPatterns
                patterns{i} = results.patternAnalysis{i}.pattern;
                probs(i) = results.patternAnalysis{i}.probability;
            end
            
            % Plot pattern probabilities
            figure;
            bar(probs);
            title('Discovered Gene Expression Patterns');
            xlabel('Pattern Index');
            ylabel('Probability');
            
            % Add pattern strings as x-tick labels
            xticks(1:numPatterns);
            xticklabels(patterns);
            xtickangle(45);
            
            % Create heatmap of top patterns
            maxPatternsToShow = min(5, numPatterns);
            
            if maxPatternsToShow > 0
                figure;
                patternMatrix = zeros(maxPatternsToShow, obj.numGenes);
                
                for i = 1:maxPatternsToShow
                    for j = 1:length(patterns{i})
                        patternMatrix(i, j) = str2double(patterns{i}(j));
                    end
                end
                
                % Use heatmap or imagesc based on MATLAB version
                try
                    h = heatmap(patternMatrix);
                    h.XDisplayLabels = obj.geneNames;
                    h.YDisplayLabels = arrayfun(@(i) ['Pattern ' num2str(i)], 1:maxPatternsToShow, 'UniformOutput', false);
                    h.Title = 'Top Gene Expression Patterns';
                    h.ColorbarVisible = 'off';
                    h.GridVisible = 'on';
                    h.ColorLimits = [0 1];
                    h.Colormap = [1 1 1; 0 0.7 0];
                catch
                    % For older MATLAB versions that don't support heatmap
                    figure;
                    imagesc(patternMatrix);
                    colormap([1 1 1; 0 0.7 0]);
                    colorbar('off');
                    axis tight;
                    title('Top Gene Expression Patterns');
                    xticks(1:obj.numGenes);
                    xticklabels(obj.geneNames);
                    yticks(1:maxPatternsToShow);
                    yticklabels(arrayfun(@(i) ['Pattern ' num2str(i)], 1:maxPatternsToShow, 'UniformOutput', false));
                    xtickangle(45);
                end
            end
        end
    end
end

% Script to demonstrate usage
% --------------------------

% Example usage
function runInverseGroverExample()
    % Generate sample expression data
    rng(42); % For reproducibility
    geneNames = {'TP53', 'MYC', 'BRCA1', 'EGFR', 'KRAS'};
    numGenes = length(geneNames);
    numSamples = 10;

    % Create random expression matrix (genes x samples)
    expressionMatrix = rand(numGenes, numSamples);

    % Create the analyzer object
    groverAnalyzer = InverseGroverGeneNetwork(expressionMatrix, geneNames);

    % Run the analysis
    [discoveredPatterns, probabilities] = groverAnalyzer.runInverseGrover();

    % Analyze the results
    results = groverAnalyzer.analyzePatterns(discoveredPatterns, probabilities);

    % Display results
    fprintf('\nAnalysis Results:\n');
    fprintf('Analyzed %d genes across %d samples\n', ...
        results.numGenesAnalyzed, results.numSamples);
    fprintf('Discovered %d significant gene expression patterns\n\n', ...
        results.numPatternsDiscovered);

    for i = 1:min(5, results.numPatternsDiscovered)
        pattern = results.patternAnalysis{i};
        fprintf('Pattern %d: %s (probability: %.4f)\n', ...
            i, pattern.pattern, pattern.probability);
        fprintf('Expressed genes: %s\n', ...
            strjoin(pattern.expressedGenes, ', '));
        fprintf('Suppressed genes: %s\n\n', ...
            strjoin(pattern.suppressedGenes, ', '));
    end

    % Visualize the results
    groverAnalyzer.visualizeResults(results);
end