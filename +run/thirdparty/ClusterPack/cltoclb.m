% function clb = cltoclb(cl)
%
% copyright (c) 1998-2011 by Alexander Strehl

function clb = cltoclb(cl)
if (prod(size(cl))==max(size(cl))),
   clb = zeros(max(cl),length(cl));
   for i=1:size(clb,1),
      clb(i,:) =  (cl==i);
   end;
else
   disp('cltoclb: cluster labels size mismatch')
end;
