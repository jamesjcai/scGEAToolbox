% function m = ceevalmutual(cls,cl),
%
% copyright (c) 1998-2011 by Alexander Strehl


function m = ceevalmutual(cls,cl),

m = 0;
totinds = 0;
n = length(cl);
for i=1:size(cls,1),
   inds = find(isfinite(cls(i,:)));
   q(i) = evalmutual(checkcl(cl(inds)),checkcl(cls(i,inds)));
   m = q(i)*length(inds) + m;
   totinds = totinds + length(inds);
end;
m = m/totinds;
