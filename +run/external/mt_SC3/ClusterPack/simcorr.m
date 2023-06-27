% function s = simcorr(a,b)
%
% DESCRIPTION
%   computes Pearson correlation based similarity between row objects
%   in matrices a and b
%
% Copyright (c) 1998-2011 by Alexander Strehl

function s = simcorr(a,b)

if ~exist('b')
  b = a;
end;

n = size(a,1);
m = size(b,1);
d = size(a,2);
if (d~=size(b,2))
  disp('simcorr: data dimensions do not match');
  return;
end;

na = a - (mean(a,2)*ones(1,size(a,2)));
nb = b - (mean(b,2)*ones(1,size(b,2)));

zaehler = na*nb';

alength = sum(na.^2,2);
blength = sum(nb.^2,2);
nenner = (alength * blength').^(1/2);

correlation = zaehler ./ nenner;
s = (correlation/2)+(1/2);

