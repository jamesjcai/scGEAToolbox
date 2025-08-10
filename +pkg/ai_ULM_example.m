
% https://www.nature.com/articles/s41525-020-00151-y
% https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformaticsadvances/2/1/10.1093_bioadv_vbac016/2/vbac016_supplementary_data.pdf?Expires=1757094748&Signature=ax6DAXDD6KbtvR9zkl6-VEZuXljEgPEXD07j69o8Gl8ya9zrvCFqFzI7AS5SZCoGERSWcaJnm7xqiUjbAKVW7JGm4n0iRkIwcm5bnowCkxQwnI4n1vW7kdRFTuWiCx7HifZmJ6o5FN5SAY~Qn76Tt7eb5pDRNFE2LcP-aq~CU3911qkue5-Ci9qDD4qx9CfVoDHU9JFPVqlEdvyiRbowt5ul3gdNWMaIgkBz2aGoijBdXK71Svcp89Pnq3mWUv0dArlz4ViWU1728yJyHdoxRSWxdSgjzFtZa1ji6lMLcgAFepTln7cltbJUckt4AAIuUX2sGqkmr9zLSPCHhuRnhg__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA

% Example usage function
function ai_ULM_example()
% Example demonstrating ULM usage

fprintf('ULM Example\n');
fprintf('===========\n');

% Create synthetic data
nGenes = 1000;
nSamples = 50;
nRegulators = 10;

% Generate expression matrix
expressionMatrix = randn(nGenes, nSamples);

% Create synthetic regulons with some biological signal
regulons = struct('name', {}, 'targets', {}, 'weights', {});

for i = 1:nRegulators
    % Random targets
    nTargets = randi([10, 30]);
    targets = randperm(nGenes, nTargets);
    weights = rand(nTargets, 1) + 0.5; % Positive weights
    
    regulons(i) = createRegulon(sprintf('TF_%d', i), targets, weights);
    
    % Add some correlated signal to targets
    activity = randn(1, nSamples);
    for j = 1:nTargets
        expressionMatrix(targets(j), :) = expressionMatrix(targets(j), :) + ...
            weights(j) * activity * 0.5;
    end
end

% Run ULM analysis
[activities, stats] = pkg.ai_ULM(expressionMatrix, regulons, 'verbose', true);

% Display results
fprintf('\nResults Summary:\n');
fprintf('Valid regulators: %d\n', stats.nValidRegulators);
fprintf('Mean activity score: %.3f ± %.3f\n', mean(stats.meanActivity), std(stats.meanActivity));

% Simple visualization
figure('Name', 'ULM Results');

subplot(2,2,1);
heatmap(activities.scores, 'Colormap', redblue, 'ColorLimits', [-2 2]);
title('Activity Scores');
xlabel('Samples');
ylabel('Regulators');

subplot(2,2,2);
histogram(activities.pvalues(:), 20);
title('P-value Distribution');
xlabel('P-value');
ylabel('Frequency');

subplot(2,2,3);
scatter(stats.meanActivity, stats.significantSamples);
title('Mean Activity vs Significant Samples');
xlabel('Mean Activity Score');
ylabel('Number of Significant Samples');

subplot(2,2,4);
bar(stats.nTargetsUsed);
title('Number of Targets per Regulator');
xlabel('Regulator Index');
ylabel('Number of Targets');

fprintf('\nExample completed. Check the figure for results.\n');

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


function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.
%   Adam Auton, 9th October 2009
if nargin < 1, m = size(get(gcf,'colormap'),1); end
if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end
c = [r g b]; 
end
