function sc_scattermarker(X,genelist,g,s)

if iscell(g)
    for k=1:length(g)
        sc_scattermarker(X,genelist,g{k},s);
    end
elseif isstring(g) && ~isStringScalar(g)
    for k=1:length(g)
        sc_scattermarker(X,genelist,g(k),s);
    end
elseif isStringScalar(g) || ischar(g)
    x=s(:,1);
    y=s(:,2);
    if ismember(g,genelist)
        z=log2(1+X(genelist==g,:));
        figure;
        sc_stemscatter(x,y,z);
        title(g)
    else
        warning('No expression.');
    end
end