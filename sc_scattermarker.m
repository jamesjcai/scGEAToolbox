function sc_scattermarker(X,genelist,g,s)

x=s(:,1);
y=s(:,2);
z=log2(1+X(genelist==g,:));
if isempty(z)
    warning('No expression.');
else
    sc_stemscatter(x,y,z);
    title(g)
end
