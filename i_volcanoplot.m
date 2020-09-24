function i_volcanoplot(fc,pvals,genelist)
x=fc;
y=-log10(pvals);
[~,idx]=maxk(abs(y),30);
scatter(x,y,5);
hold on
% scatter(x(idx),y(idx),'rx');
for k=1:length(idx)
    text(x(idx(k))+0.05,y(idx(k)),genelist(idx(k)));
end
