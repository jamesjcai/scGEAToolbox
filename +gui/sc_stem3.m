function sc_stem3(X, Y, genelist, numgene)

if nargin < 4, numgene = 50; end
if ~isempty(Y), assert(size(X, 1) == size(Y, 1)); end

nz = size(X, 2);
bigX = zeros(numgene, nz);
bigY = zeros(numgene, nz);
bigZ = zeros(numgene, nz);
for k = 1:numgene
    bigX(k, :) = k;
    bigY(k, :) = 1:nz;
    bigZ(k, :) = X(k, :);
end

stem3(bigX, bigY, bigZ, 'marker', 'none');
if ~isempty(Y)
    ny = size(Y, 2);
    bigX = zeros(numgene, ny);
    bigY = zeros(numgene, ny);
    bigZ = zeros(numgene, ny);
    for k = 1:numgene
        bigX(k, :) = k;
        bigY(k, :) = nz + (1:ny);
        bigZ(k, :) = Y(k, :);
    end
    hold on
    stem3(bigX, bigY, bigZ, 'marker', 'none');
end
view(60, 35)
xlabel('Genes');
ylabel('Cells');
zlabel('Expression');
set(gca, 'XTick', 1:numgene);
set(gca, 'XTickLabel', genelist(1:numgene));


