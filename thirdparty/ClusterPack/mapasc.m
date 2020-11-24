% function din = mapasc(din,column)
%
% copyright (c) 1998-2011 by Alexander Strehl

function din = mapasc(din,column)
if ~exist('column')
  column = 1;
end;
last = NaN;
curind = 0;
for i=1:size(din,1),
  if (last~=din(i,column))|(isnan(last)),
    last = din(i,column);
    curind = curind + 1;
  end;
  din(i,column) = curind;
end;
