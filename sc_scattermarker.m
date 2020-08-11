function sc_scattermarker(X,genelist,g,s,methodid)
%SC_SCATTERMARKER(X,genelist,g,s,methodid)
%
% USAGE:
% s=sc_tsne(X,3);
% g=["AGER","SFTPC","SCGB3A2","TPPP3"];
% sc_scattermarker(X,genelist,g,s);


if nargin<5, methodid=1; end
if iscell(g)
    for k=1:length(g)
        sc_scattermarker(X,genelist,g{k},s,methodid);
    end
elseif isstring(g) && ~isStringScalar(g)
    for k=1:length(g)
        sc_scattermarker(X,genelist,g(k),s,methodid);
    end
elseif isStringScalar(g) || ischar(g)
    if ismember(g,genelist)
        x=s(:,1);
        y=s(:,2);            
        figure;        
        switch methodid
            case 1
                z=log2(1+X(genelist==g,:));
                sc_stemscatter(x,y,z);
            case 2
                z=s(:,3);
                c=log2(1+X(genelist==g,:));
                scatter3(x,y,z,10,c,'filled');
        end
        title(g)
    else
        warning('%s no expression',g);
    end
end
