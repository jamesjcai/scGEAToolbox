function [X, genelist] = sc_rmdugenes(X, genelist)

    % Remove genes with duplicate name - Optimized version
    [~, first_idx] = unique(genelist, 'stable');

    if length(first_idx) == length(genelist)
        return;
    else
        num_duplicates = length(genelist) - length(first_idx);
        warning('Duplicate gene names are removed. %d duplicate genes found.', num_duplicates);
        X = X(first_idx, :);
        genelist = genelist(first_idx);
    end

end

%{
 Pre-allocate result matrix
X_new = zeros(length(first_idx), size(X, 2));

% Vectorized summation for each unique gene group
for i = 1:length(first_idx)
    group_members = (first_idx == i);
    X_new(i, :) = sum(X(group_members, :), 1);
end

% Update outputs
X = X_new;
genelist = genelist(first_idx);
%}



%{
function [X, genelist] = sc_rmdugenes(X, genelist)
% Remove genes with duplicate name

[~, w] = unique(genelist, 'stable');
duplicate_indices = setdiff(1:numel(genelist), w);
if isempty(duplicate_indices)
    return;
end

for k = 1:length(duplicate_indices)
    idx = find(genelist == genelist(duplicate_indices(k)));
    X(idx(1), :) = sum(X(idx, :), 1);
end
X(duplicate_indices, :) = [];
genelist(duplicate_indices) = [];
end
%}

%{
function [X, genelist] = sc_rmdugenes(X, genelist)
% Remove genes with duplicate name - Optimized version
[~, first_idx, group_idx] = unique(genelist, 'stable');

% If no duplicates, return early
if length(first_idx) == length(genelist)
    return;
end

% Pre-allocate result matrix
X_new = zeros(length(first_idx), size(X, 2));

% Vectorized summation for each unique gene group
for i = 1:length(first_idx)
    group_members = (group_idx == i);
    if sum(group_members) == 1
        X_new(i, :) = X(group_members, :);
    else
        X_new(i, :) = sum(X(group_members, :), 1);
    end
end

% Update outputs
X = X_new;
genelist = genelist(first_idx);
end
%}