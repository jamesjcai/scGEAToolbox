function [Tct] = i_celltypebrushed(X, genelist, s, ...
    brushedData, species, ~, database, bestonly, ~)

if nargin < 9, subtype = 'all'; end
if nargin < 8, bestonly = true; end
if nargin < 7, database = 'panglaodb'; end
if nargin < 6, organ = 'all'; end
if nargin < 5, species = 'mouse'; end

if islogical(brushedData)
    i = brushedData;
else
    [~, i] = ismember(brushedData, s, 'rows');
end
Xi = X(:, i);
gi = upper(genelist);
if strcmpi(database, 'clustermole')
    [Tct] = run.r_clustermole(Xi, gi, [], 'species', species);
elseif strcmpi(database, 'panglaodb')


    [Tct] = run.ml_alona(Xi, gi, [], 'species', species, ...
        'bestonly', bestonly);


end
end
