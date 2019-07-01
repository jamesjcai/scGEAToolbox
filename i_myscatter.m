function h=i_myscatter(s,c)
x=s(:,1);
y=s(:,2);

if nargin<2
    c='b';
end

if size(s,2)==3
    z=s(:,3);
end

if nargout>0
    if size(s,2)==2
                h=scatter(x,y,10,c,'filled');
    elseif size(s,2)==3
                h=scatter3(x,y,z,10,c,'filled');
    end
else
    if size(s,2)==2
                scatter(x,y,10,c,'filled');
    elseif size(s,2)==3
                scatter3(x,y,z,10,c,'filled');
    end   
end
