nchoosek(20, 3)
X = zeros(35);
%%
[X,Y] = meshgrid(1:35, 1:35);
Z = reshape(randn(35.^2, 1), 35, 35);
surface(X,Y,Z)
view(3)

%%
C = nchoosek(1:20, 3);
v = randn(length(C), 1);

numPoints = length(C);
dataset = C; % Random integers 0, 1, or 2

customDistance = @(x, y) 3 - sum(x == y, 2);

distanceMatrix = zeros(numPoints); % Preallocate matrix
for i = 1:numPoints
    distanceMatrix(i, :) = arrayfun(@(j) customDistance(dataset(i, :), dataset(j, :)), 1:numPoints);
end

% Step 4: Build KNN graph
k = 6; % Number of nearest neighbors
knnGraph = zeros(numPoints); % Adjacency matrix for KNN graph
for i = 1:numPoints
    [~, sortedIndices] = sort(distanceMatrix(i, :)); % Sort by distance
    nearestNeighbors = sortedIndices(2:k+1); % Exclude itself (distance=0)
    knnGraph(i, nearestNeighbors) = 1; % Mark edges in adjacency matrix
end


g=digraph(knnGraph, string(1:1140)');
gui.i_singlegraph(g)