% function s = simcosi(a,b)
%
% DESCRIPTION
%   computes cosine based similarity between row objects in matrices a and b
%
% Copyright (c) 1998-2011 by Alexander Strehl

function s = simcosi(a,b)

if ~exist('b')
  b = a;
end;

n = size(a,1);
m = size(b,1);
d = size(a,2);
if (d~=size(b,2))
  disp('simcosi: data dimensions do not match');
  return;
end;

alength = (sum((a.^2),2)).^(1/2); % one column
blength = (sum((b.^2),2)).^(1/2); % one column
s = (a * b') ./ ( alength * blength' );


