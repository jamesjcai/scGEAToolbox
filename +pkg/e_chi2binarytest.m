function [chi2_stat, p_value, contingency_table, expected_freq] = e_chi2binarytest(var1, var2, alpha)
% CHI2_BINARY_TEST Performs chi-squared test of independence for two binary variables
%
% Inputs:
%   var1 - binary vector (0s and 1s) for first variable
%   var2 - binary vector (0s and 1s) for second variable
%   alpha - significance level (default: 0.05)
%
% Outputs:
%   chi2_stat - chi-squared test statistic
%   p_value - p-value of the test
%   contingency_table - 2x2 contingency table
%   expected_freq - expected frequencies under null hypothesis

if nargin < 3
    alpha = 0.05;
end

% Ensure inputs are column vectors
var1 = var1(:);
var2 = var2(:);

% Check that variables are binary
if ~all(ismember(var1, [0, 1])) || ~all(ismember(var2, [0, 1]))
    error('Both variables must be binary (contain only 0s and 1s)');
end

% Check equal lengths
if length(var1) ~= length(var2)
    error('Variables must have the same length');
end

% Count the four cells directly rather than calling crosstab and padding
% its result. crosstab only returns rows and columns for levels that
% actually occur, so a constant variable gave a 1-by-2 (or 2-by-1) table
% and the padding appended the zero row or column at the END. When the
% constant value was 1, every count therefore landed in the var1=0 row (or
% the var2=0 column) and the printed table said the opposite of the truth:
% for 200 cells with var2 all ones, it showed 200 in the var2=0 column and
% 0 in var2=1.
contingency_table = [ ...
    sum(var1 == 0 & var2 == 0), sum(var1 == 0 & var2 == 1); ...
    sum(var1 == 1 & var2 == 0), sum(var1 == 1 & var2 == 1)];

% Calculate expected frequencies
row_totals = sum(contingency_table, 2);
col_totals = sum(contingency_table, 1);
n = sum(contingency_table(:));

expected_freq = (row_totals * col_totals) / n;

% Degrees of freedom for 2x2 table
df = 1;

% A constant variable leaves a whole row or column of the table empty, so
% the matching expected frequencies are zero and the statistic is 0/0.
% That used to reach the verdict block below as NaN, and `if NaN < alpha`
% is false, so the function printed "FAIL TO REJECT" and "No significant
% association between variables" -- a confident negative conclusion drawn
% from an undefined statistic. That text is shown to the user verbatim by
% +gui/callback_2GeneCooccurrenceTest, which captures it with evalc, and
% asking about two genes where one is simply not detected in the selected
% cells is an ordinary thing to do.
isDegenerate = any(row_totals == 0) || any(col_totals == 0);
if isDegenerate
    chi2_stat = NaN;
    p_value = NaN;
    warning('pkg:e_chi2binarytest:constantVariable', ...
        ['At least one variable takes a single value, so the ', ...
        'contingency table has an empty row or column and the ', ...
        'chi-squared statistic is undefined. No test was performed.']);
else
    % Calculate chi-squared statistic
    chi2_stat = sum(sum((contingency_table - expected_freq).^2 ./ expected_freq));

    % Calculate p-value
    p_value = 1 - chi2cdf(chi2_stat, df);
end

% Display results
fprintf('\n=== Chi-squared Test for Independence ===\n');
fprintf('Null hypothesis: Variables are independent\n');
fprintf('Alternative hypothesis: Variables are dependent\n\n');

fprintf('Contingency Table:\n');
fprintf('         var2=0  var2=1  Total\n');
fprintf('var1=0   %6d  %6d  %6d\n', contingency_table(1,1), contingency_table(1,2), sum(contingency_table(1,:)));
fprintf('var1=1   %6d  %6d  %6d\n', contingency_table(2,1), contingency_table(2,2), sum(contingency_table(2,:)));
fprintf('Total    %6d  %6d  %6d\n', sum(contingency_table(:,1)), sum(contingency_table(:,2)), n);

fprintf('\nExpected Frequencies:\n');
fprintf('         var2=0  var2=1\n');
fprintf('var1=0   %6.2f  %6.2f\n', expected_freq(1,1), expected_freq(1,2));
fprintf('var1=1   %6.2f  %6.2f\n', expected_freq(2,1), expected_freq(2,2));

fprintf('\nTest Statistics:\n');
fprintf('Chi-squared statistic: %.4f\n', chi2_stat);
fprintf('Degrees of freedom: %d\n', df);
fprintf('P-value: %.6f\n', p_value);
fprintf('Significance level: %.3f\n', alpha);

if isnan(p_value)
    % Never fall through to a verdict on a NaN.
    fprintf('Result: NO TEST PERFORMED\n');
    fprintf(['Conclusion: at least one variable is constant, so ', ...
        'independence cannot be assessed.\n']);
elseif p_value < alpha
    fprintf('Result: REJECT null hypothesis (p < %.3f)\n', alpha);
    fprintf('Conclusion: Variables are significantly associated\n');
else
    fprintf('Result: FAIL TO REJECT null hypothesis (p >= %.3f)\n', alpha);
    fprintf('Conclusion: No significant association between variables\n');
end

% Check assumptions
min_expected = min(expected_freq(:));
if ~isDegenerate && min_expected < 5
    fprintf('\nWarning: Minimum expected frequency (%.2f) is less than 5.\n', min_expected);
    fprintf('Chi-squared test may not be appropriate. Consider Fisher''s exact test.\n');
end

end

% Example usage and demonstration
function demo_chi2_binary()
fprintf('=== DEMONSTRATION ===\n');

% Example 1: Independent variables
fprintf('\nExample 1: Testing independent variables\n');
% Save and restore the caller's random stream. Seeding this demonstration's data is
% fine; leaving the session parked on that seed is not -- it then
% governs every later tsne, umap and clustering call in the session.
rngState = rng();
restoreRng = onCleanup(@() rng(rngState));
rng(42); % For reproducibility
n = 200;
var1 = binornd(1, 0.3, n, 1);  % 30% probability of 1
var2 = binornd(1, 0.4, n, 1);  % 40% probability of 1, independent

[chi2_stat1, p_val1] = pkg.e_chi2binarytest(var1, var2);

% Example 2: Dependent variables
fprintf('\n\nExample 2: Testing dependent variables\n');
var1 = binornd(1, 0.3, n, 1);
% Make var2 dependent on var1
var2 = zeros(n, 1);
var2(var1 == 1) = binornd(1, 0.8, sum(var1), 1);  % Higher prob when var1=1
var2(var1 == 0) = binornd(1, 0.2, sum(var1==0), 1);  % Lower prob when var1=0

[chi2_stat2, p_val2] = pkg.e_chi2binarytest(var1, var2);

fprintf('\n=== SUMMARY ===\n');
fprintf('Example 1 (independent): chi2=%.4f, p=%.6f\n', chi2_stat1, p_val1);
fprintf('Example 2 (dependent):   chi2=%.4f, p=%.6f\n', chi2_stat2, p_val2);
end

% Run demonstration
% demo_chi2_binary();
