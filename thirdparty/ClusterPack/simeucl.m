% function s = simeucl(a,b)
%
% DESCRIPTION
%   computes Euclidean based similarity between row objects in matrices a and b
%
% Copyright (c) 1998-2011 by Alexander Strehl

function s = simeucl(a,b)

if ~exist('b')
  b = a;
end;

n = size(a,1);
m = size(b,1);
d = size(a,2);
if (d~=size(b,2))
  disp('simeucl: data dimensions do not match');
  return;
end;        

s = exp(-((ones(m, 1) * sum((a.^2)', 1))' + ones(n, 1) * sum((b.^2)',1) - 2.*(a*(b'))));
