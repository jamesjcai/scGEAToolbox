% function s = simbjac(a,b)
%
% DESCRIPTION
%   computes Jaccard similarity between row objects in matrices a and b
%
% Copyright (c) 1998-2011 by Alexander Strehl

function s = simbjac(a,b)

if ~exist('b')
  b = a;
end;

s = simxjac(a>0,b>0);

