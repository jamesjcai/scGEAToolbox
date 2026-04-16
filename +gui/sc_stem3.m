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
% figure;

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


%                 figure;
%                 stem3(app.X(idx(1:50),:),'marker','none');
%                 view(139,35)
%                 ylabel('Genes');
%                 xlabel('Cells');
%                 zlabel('Expression');
%                 set(gca,'YTick',1:30);
%                 set(gca,'YTickLabel',table2array(app.T(1:30,1)));
%                 ylim([0 31]);
