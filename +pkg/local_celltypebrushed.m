function [Tct] = local_celltypebrushed(X, genelist, s, ...
    brushedData, species, organ, database, bestonly, subtype)

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
%[Xi,gi]=sc_selectg(Xi,genelist);
if strcmpi(database, 'clustermole')
    %disp('Using clustermole marker database')
    [Tct] = run.r_clustermole(Xi, gi, [], 'species', species);
elseif strcmpi(database, 'panglaodb')
    % disp('Using panglaodb marker database')

    
%    [Tct] = run.mt_alona(Xi, gi, [], 'species', species, 'organ', organ, ...
%        'bestonly', bestonly, 'subtype', subtype);
    
    
    [Tct] = run.mt_alona_new(Xi, gi, [], 'species', species, ...
        'bestonly', bestonly);
    

end
end
