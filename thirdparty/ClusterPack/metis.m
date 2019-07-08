% function labels=metis(x,k)
%
% copyright (c) 1998-2011 by Alexander Strehl

function labels=metis(x,k) 

filename = wgraph(x,[],0);
labels = sgraph(k,filename);
delete(filename);
