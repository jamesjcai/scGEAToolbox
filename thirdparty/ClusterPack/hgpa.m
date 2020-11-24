% function cl = hgpa(cls,k)
%
% DESCRIPTION
%  Performs HGPA for CLUSTER ENSEMBLES
%
% Copyright (c) 1998-2011 by Alexander Strehl

function cl = hgpa(cls,k)

disp('CLUSTER ENSEMBLES using HGPA');

if ~exist('k'),
   k = max(max(cls));
end;

r = size(cls,1);
clb = [];
for i=1:r,
   clb = [clb; cltoclb(cls(i,:))];
   kq(i) = max(cls(i,:));
   lastindex = sum(kq(1:i));
   firstindex = lastindex-kq(i)+1;
   xy(firstindex:lastindex,:) =  [i*ones(kq(i),1) ((1:kq(i))')];
end;

cl = clhgraph(clb',k,ones(1,size(clb,1))); 
