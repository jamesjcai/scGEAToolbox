% function ai_AUCell_example()
% Example of how to use the AUCell function

fprintf('AUCell Example\n');
fprintf('==============\n');

% Create example data
nGenes = 2000;
nCells = 100;
expressionMatrix = rand(nGenes, nCells) .* 10;

% Add some signal to specific genes for demonstration
signalGenes1 = 1:50;  % First gene set
signalGenes2 = 101:180;  % Second gene set

% Make some cells have higher expression of gene set 1
expressionMatrix(signalGenes1, 1:30) = expressionMatrix(signalGenes1, 1:30) + 5;

% Make some cells have higher expression of gene set 2
expressionMatrix(signalGenes2, 31:60) = expressionMatrix(signalGenes2, 31:60) + 5;

% Define gene sets
geneSets = {signalGenes1, signalGenes2, 201:250};  % Third set as control

% Run AUCell
[aucScores, rankings] = pkg.ai_AUCell(expressionMatrix, geneSets);

% Display results
fprintf('\nAUC Scores Summary:\n');
for i = 1:size(aucScores, 1)
    fprintf('Gene Set %d: mean=%.3f, std=%.3f\n', i, mean(aucScores(i,:)), std(aucScores(i,:)));
end

% Simple visualization of results
figure('Name', 'AUCell Results');
for i = 1:size(aucScores, 1)
    subplot(size(aucScores, 1), 1, i);
    bar(aucScores(i, :));
    title(sprintf('Gene Set %d AUC Scores', i));
    xlabel('Cells');
    ylabel('AUC Score');
    ylim([0, 1]);
end

fprintf('\nExample completed. Check the figure for results.\n');

%end