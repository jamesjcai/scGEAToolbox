function h=i_myscatter(s,c,methodid)

if nargin<3, methodid=1; end
if nargin<2, c='b'; end

x=s(:,1);
y=s(:,2);

if iscell(c), c=grp2idx(c); end
isbinary=logical(numel(unique(c))==2);

switch methodid
    case 1
        if size(s,2)==3, z=s(:,3); end
        if size(s,2)==2
           h=scatter(x,y,5,c,'filled');
        elseif size(s,2)==3
           h=scatter3(x,y,z,5,c,'filled');
        end
    case 2
        h=gscatter(x,y,c,[],[],10);
end

if isbinary
    colormap(lines(2));
end
