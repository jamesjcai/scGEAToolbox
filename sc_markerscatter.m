function sc_markerscatter(X,genelist,g,s,methodid,sz)
%SC_MARKERSCATTER(X,genelist,g,s,methodid)
%
% USAGE:
% s=sc_tsne(X,3);
% g=["AGER","SFTPC","SCGB3A2","TPPP3"];
% sc_scattermarker(X,genelist,g,s);

if isvector(s)||isscalar(s), error('S should be a matrix.'); end
if nargin<6, sz=5; end
if nargin<5, methodid=1; end
if iscell(g)
    for k=1:length(g)
        figure;
        sc_scattermarker(X,genelist,g{k},s,methodid,sz);
    end
elseif isstring(g) && ~isStringScalar(g)
    for k=1:length(g)
        figure;
        sc_scattermarker(X,genelist,g(k),s,methodid,sz);
    end
elseif isStringScalar(g) || ischar(g)
    if ismember(g,genelist)
        x=s(:,1);
        y=s(:,2);
        if min(size(s))==2
            z=[];
        else
            z=s(:,3);
        end               
        switch methodid
            case 1
                z=log2(1+X(genelist==g,:));
                sc_stemscatter(x,y,z);
            case 2                
                c=log2(1+X(genelist==g,:));
                if isempty(z)
                    scatter(x,y,sz,c,'filled');
                else
                    scatter3(x,y,z,sz,c,'filled');
                end
                colormap('default');
            case 3                
                c=log2(1+X(genelist==g,:));
                if isempty(z)
                    scatter(x,y,sz,c,'filled');
                else
                    scatter3(x,y,z,sz,c,'filled');
                end                
                a=colormap('autumn');
                a(1,:)=[.8 .8 .8];
                colormap(a);
        end
        title(sprintf('%s %f%%',g,100*sum(c>0)./numel(c)));
    else
        warning('%s no expression',g);
    end
end
