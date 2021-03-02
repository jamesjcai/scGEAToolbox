% function cl = clagglmin(x,k,sfct)
%
% DESCRIPTION
%   provides cluster labels 1 to k using single link agglomerative clustering
%   this is a slow O(n^3) implementation
%
% Copyright (c) 1998-2011 by Alexander Strehl


function cl = clagglmin(x,k,sfct)

if ~exist('sfct'),
   sfct = 'simeucl';
   disp('clagglmin: using default Euclidean measure');
end;

n = size(x,1);
s = feval(sfct,x);
disp('clagglmin: computed s');
for i=1:size(s,1),
  s(i,i) = -inf;
end;

cl = 1:n;
for i=1:(n-k),
  if ((n-k>500)&(mod(n-i+1,100)==0)), disp(['clagglmin: at ' num2str(n-i+1) ' clusters ']); end;
  % compute indices of closest pair
  %[a,b] = ind2sub(size(s),topnindx(s,1)); % less efficient
  [mm, ii] = max(reshape(s,1,prod(size(s))));
  [a,b] = ind2sub(size(s),ii);
  % merge clusters by replacing all members of a with b's label
  MembersOfAInd = find(cl==(cl(a)));
  MembersOfBInd = find(cl==(cl(b)));
  cl(MembersOfAInd) = ones(1,length(MembersOfAInd)) * cl(b);
  % disable pairwise similarity
  blankblock = -inf * ones(length(MembersOfAInd),length(MembersOfBInd));
  s(MembersOfAInd,MembersOfBInd) = blankblock;
  s(MembersOfBInd,MembersOfAInd) = blankblock';
end;
cl = onetomax(cl);

