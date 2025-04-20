    % Generate sample expression data
    % rng(42); % For reproducibility
    geneNames = {'TP53', 'MYC', 'BRCA1', 'EGFR', 'KRAS'};
    numGenes = length(geneNames);
    numSamples = 10;

    % Create random expression matrix (genes x samples)
    expressionMatrix = rand(numGenes, numSamples);

    % Create the analyzer object
    groverAnalyzer = qtm.InverseGroverGeneNetwork(expressionMatrix, geneNames);

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