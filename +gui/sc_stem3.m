function sc_stem3(X, Y, genelist, numgene)

if nargin < 4, numgene = 50; end
if ~isempty(Y), assert(size(X, 1) == size(Y, 1)); end

bigX = [];
bigY = [];
bigZ = [];
for k = 1:numgene
    z = X(k, :);
    x = k * ones(size(z));
    y = 1:length(z);
    bigX = [bigX; x];
    bigY = [bigY; y];
    bigZ = [bigZ; z];
end
nz = length(z);
% figure;

stem3(bigX, bigY, bigZ, 'marker', 'none');
if ~isempty(Y)
    bigX = [];
    bigY = [];
    bigZ = [];
    for k = 1:numgene
        z = Y(k, :);
        x = k * ones(size(z));
        y = nz + (1:length(z));
        bigX = [bigX; x];
        bigY = [bigY; y];
        bigZ = [bigZ; z];
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
