function [h]=i_gscatter3(s,c,methodid,targetc)

if nargin<4, targetc=1; end
if nargin<3, methodid=1; end
if nargin<2, c=ones(size(s,1),1); end

x=s(:,1);
y=s(:,2);

if iscell(c)||isstring(c)
    c=grp2idx(c);
end
kc=numel(unique(c));

if size(s,2)>=3, z=s(:,3); end

switch methodid
    case 1
        if size(s,2)==2
           h=scatter(x,y,10,c);
           % box on
        elseif size(s,2)>=3
           h=scatter3(x,y,z,10,c);
%            i=c==1;
%            h=scatter3(x(i),y(i),z(i),10);
%            hold on
%            for k=2:max(c)
%                i=c==k;
%                scatter3(x(i),y(i),z(i),10);
%            end
        end
    case 2
        if size(s,2)==2
            h=gscatter(x,y,c,[],[],5);
        elseif size(s,2)>=3
            h=gscatter3b(x,y,z,c,[],[],5);
        end
        box off
    case 3
        h=gui.i_gscatter3(s,c,1);
        h.MarkerEdgeAlpha=0;
        hold on
        i=c==targetc;
        h=gui.i_gscatter3(s(i,:),c(i),1);
        hold off        
end

if kc<=50
    colormap(lines(kc));
else
%     cx=colormap('autumn');
%     cx(1,:)=[.8 .8 .8];
%     colormap(cx);
    colormap default
end
grid on

end


