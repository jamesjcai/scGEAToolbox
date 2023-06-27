% function mapout = fastchangem(map,newcode,oldcode)
%
% copyright (c) 1998-2011 by Alexander Strehl

function mapout = fastchangem(map,newcode,oldcode)

mapout=map;
for i=1:length(newcode),
  mapout(find(map==oldcode(i))) = newcode(i);
end;
