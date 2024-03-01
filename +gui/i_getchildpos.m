function [px_new, childfigpos] = i_getchildpos(partentfigpos, childfigpos)

px_new=[];
if ismac, return; end

if ~isequal(size(partentfigpos),[1 4])
    partentfigpos = partentfigpos.Position;
end
if ~isequal(size(childfigpos),[1 4])
    childfigpos = childfigpos.Position;
end

    p = partentfigpos;
    cx = [p(1)+p(3)/2 p(2)+p(4)/2];

    px = childfigpos;
    px_new = [cx(1)-px(3)/2 cx(2)-px(4)/2];
if nargout>1
    childfigpos([1 2]) = px_new; 
end