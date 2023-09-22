function i_dotplot_old(X0, X1, genelist, tgene, uselog)

if nargin < 5
    uselog = false;
end

[i] = ismember(tgene, genelist);
if ~any(i), return; end
z = length(tgene) - sum(i);
if z > 0, fprintf('%d gene(s) not in the list are excluded.\n', z); end
tgene = tgene(i);

% tgene=string(T.gene(1:10));
y = (1:length(tgene))';
x = [-ones(size(y)); ones(size(y))] ./ 2;
y = [y; y];


sz = ones(size(y));
c = ones(size(y));
for k = 1:length(tgene)
    a0 = X0(genelist == tgene(k), :);
    a1 = X1(genelist == tgene(k), :);
    sz(k) = sum(a0 ~= 0) ./ length(a0);
    sz(k+length(tgene)-1) = sum(a1 ~= 0) ./ length(a1);
    c(k) = mean(a0);
    c(k+length(tgene)-1) = mean(a1);
end

%c=linspace(1,10,length(x));
%c=zscore(c);
%c=c+min(c);
if uselog
    c = log2(1+c);
end
txgene = [" "; tgene];

% figure;
%sz=randi(100,1,length(x));
%scatter([-.5 .5],[-1 -1],[1 500],'k','filled');
%hold on
% x;
% y;
% sz;
% c;
scatter(x, y, 500*sz, c, 'filled');
hold on
scatter(x, y, 500*sz, 'k');
xlim([-1, 1]);
ylim([0, length(txgene)]);
colorbar
set(gca, 'YTick', 0:length(tgene))
set(gca, 'YTickLabel', txgene)
set(gca, 'XTickLabel', {'', 'WT', '', 'KO', ''})
colormap(flipud(bone));
box on
grid on
hFig = gcf;
hFig.Position(3) = hFig.Position(3) * 0.7;
