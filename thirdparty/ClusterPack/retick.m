% function retick(h,xtxt,ytxt),
%
% copyright (c) 1998-2011 by Alexander Strehl

function retick(h,xtxt,ytxt),

if ~isempty(xtxt) 
   set(h,'XTickMode','manual');
   set(h,'XTickLabelMode','manual');
   set(h,'XTick',1:size(xtxt,1));
   set(h,'XTickLabel',xtxt);
end;

if ~isempty(ytxt)
   set(h,'YTickMode','manual');
   set(h,'YTickLabelMode','manual');
   set(h,'YTick',1:size(ytxt,1));
   set(h,'YTickLabel',ytxt);
end;

axis off;
axis on;
