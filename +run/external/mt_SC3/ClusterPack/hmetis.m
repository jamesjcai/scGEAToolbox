% function labels=hmetis(x,k,w)
%
% copyright (c) 1998-2011 by Alexander Strehl

function labels=hmetis(x,k,w) 

if ~exist('w'),
  filename = wgraph(x,[],2);
else
  filename = wgraph(x,w,3);
end; 
labels = sgraph(k,filename); 
delete(filename);
