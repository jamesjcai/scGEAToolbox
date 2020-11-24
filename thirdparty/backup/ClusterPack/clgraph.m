% function cl = clgraph(x,k,sfct)
% 
% DESCRIPTION
%   provides cluster labels 1 to k from edge weighted graph partitioning 
%
% Copyright (c) 1998-2011 by Alexander Strehl


function cl = clgraph(x,k,sfct)

cl = metis(checks(feval(sfct,x)),k);
