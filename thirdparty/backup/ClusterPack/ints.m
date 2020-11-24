% function s = ints(s,tafter)
%
% copyright (c) 1998-2011 by Alexander Strehl

function s = ints(s,tafter)

if ~exist('tafter'),
  tafter = 100000000;
end;

tbefore = sum(sum(s));
if (tbefore~=0),
  s = round(s.*(tafter/tbefore));
end;

p = 1-sum(sum(s>0))/prod(size(s)); 


