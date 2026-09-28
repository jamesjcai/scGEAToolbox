function [T, Tnt] = in_DVTableProcess(T, ~, ~, direction, numPermutations)
% DIRECTION is the SC_DVG option that produced DiffSign ('mean', the
% default, or 'deviation'); it only changes how the Note describes DiffSign.
% NUMPERMUTATIONS is SC_DVG's NumPermutations (default 0); it only changes
% how the Note describes pval.

if nargin < 4 || isempty(direction), direction = 'mean'; end
if nargin < 5 || isempty(numPermutations), numPermutations = 0; end

if nargout>1
    switch direction
        case 'deviation'
            signText = ['Sign of difference in deviation from curve ' ...
                '(+1: more variable in sample 1)'];
        otherwise
            signText = ['Sign of difference in mean expression ' ...
                '(+1: up-regulated, higher mean in sample 1)'];
    end
    if numPermutations > 0
        pvalText = sprintf(['p-value of DV test (permutation null, %d ' ...
            'relabellings of the cells)'], numPermutations);
    else
        pvalText = 'p-value of DV test';
    end

    Item = T.Properties.VariableNames';

    Description = {'gene name';'log mean in sample 1';...
        'log CV in sample 1'; 'dropout rate in sample 1';...
        'distance to curve 1';'p-value of distance in sample 1';...
        'FDR of distance in sample 1';'log mean in sample 2';...
        'log CV in sample 2'; 'dropout rate in sample 2';...
        'distance to curve 2'; 'p-value of distance in sample 2';...
        'FDR of distance in sample 2'; 'Difference in distances';...
        signText; pvalText};

    if length(Item) == length(Description)
        Tnt = table(Item, Description);
    else
        assignin("base","Item", Item);
        assignin("base","Description", Description);
        Tnt = table(Item);
        warning('Variables must have the same number of rows.');
    end
end

end
