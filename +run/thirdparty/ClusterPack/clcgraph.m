% function cl = clcgraph(x,k,sfct)
%
% DESCRIPTION
%   provides cluster labels 1 to k from edge weighted value balanced
%   graph partitioning
% BUGS 
%   pmetis doesn't seem to work on PCWIN for vertex weighted graphs
%
% Copyright (c) 1998-2011 by Alexander Strehl


function cl = clcgraph(x,k,sfct)

cl = cmetis(checks(feval(sfct,x)),sum(x,2),k);
