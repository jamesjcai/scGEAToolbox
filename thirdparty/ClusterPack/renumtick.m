% function renumtick(h,xcl,ycl)
%
% copyright (c) 1998-2011 by Alexander Strehl

function renumtick(h,xcl,ycl)

if exist('xcl')&~isempty(xcl),
   xtxt = num2str(perminv(xcl)');
else
   xtxt = [];
end;
if exist('ycl')&~isempty(ycl),
   ytxt = num2str(perminv(ycl)');
else
   ytxt = [];
end;
retick(h,xtxt,ytxt);

