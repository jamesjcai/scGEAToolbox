% function cl = clhgraph(x,k,sfct)
%
% DESCRIPTION
%   provides cluster labels 1 to k from hypergraph partitioning 
%   sfct is ignored
%
% Copyright (c) 1998-2011 by Alexander Strehl


function cl = clhgraph(x,k,sfct)

cl = hmetis(x,k);
