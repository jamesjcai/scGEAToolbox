% function x = clucent(cpm,cl)
%
% Copyright (c) 1998-2011 by Alexander Strehl


function x = clucent(cpm,cl)

x = sparse(max(cl),size(cpm,2));
for i=1:max(cl)
  x(i,:) = mean(cpm(find(cl==i),:),1);
end;

