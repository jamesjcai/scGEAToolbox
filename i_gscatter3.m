function [h]=i_gscatter3(s,c,methodid)

if nargin<3, methodid=1; end
if nargin<2, c=ones(size(s,1),1); end

x=s(:,1);
y=s(:,2);

if iscell(c), c=grp2idx(c); end
kc=numel(unique(c));

if size(s,2)>=3, z=s(:,3); end

switch methodid
    case 1
        if size(s,2)==2
           h=scatter(x,y,10,c);
        elseif size(s,2)>=3
           h=scatter3(x,y,z,10,c);
        end
    case 2
        if size(s,2)==2
            h=gscatter(x,y,c,[],[],5);
        elseif size(s,2)>=3
            h=gscatter3b(x,y,z,c,[],[],5);
        end
end

if kc<=50
    colormap(lines(kc));
else
    cx=colormap('autumn');
    cx(1,:)=[.8 .8 .8];
    colormap(cx);
end
% add_3dcamera;
end


