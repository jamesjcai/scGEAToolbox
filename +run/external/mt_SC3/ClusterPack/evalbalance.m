% function p = evalbalance(trueclass,cl,x,sfct)
% 
% DESCRIPTION
%   computes balance of clustering cl
%   ignores trueclass, x, and sfct - just for compatibility
%
% Copyright (c) 1998-2011 by Alexander Strehl

function p = evalbalance(trueclass,cl,x,sfct)

o = hist(cl,1:max(cl));
imbalance = length(o)./length(cl)*max(o);
p = 1/imbalance;
