% function labels=cmetis(e,w,k)
%
% copyright (c) 1998-2011 by Alexander Strehl


function labels=cmetis(e,w,k) 

filename = wgraph(e,w,1);
labels = sgraph(k,filename);
delete(filename);
