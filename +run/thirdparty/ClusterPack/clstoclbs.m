% function clbs = clstoclbs(cls)
%
% copyright (c) 1998-2011 by Alexander Strehl

function clbs = clstoclbs(cls)

clbs = [];
for i=1:size(cls,1),
  clbs = [clbs; cltoclb(cls(i,:))];
end;
