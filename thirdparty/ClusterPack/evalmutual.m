% function p = evalmutual(trueclass,cl,x,sfct)
%
% DESCRIPTION
%   computes mutual information (transinformation) between
%   category labels trueclass and cluster labels cl
%   normalized by geometric mean of entropies
%   ignores x and sfct - just for compatibility
%
% Copyright (c) 1998-2011 by Alexander Strehl

function p = evalmutual(trueclass,cl,x,sfct)

remappedcl = zeros(size(cl));
A = zeros(max(cl),2+max(trueclass));
for i=1:max(cl),
  activepoints = find(cl==i);
  composition = hist(trueclass(activepoints),1:max(trueclass));
  j = find(composition==max(composition));
  j = j(1);
  A(i,:) = [j i composition]; 
end;
A = sortrows(A);

A = A(:,3:size(A,2));

pdf = A./sum(sum(A));
px = sum(pdf,1);
py = sum(pdf,2);

if length(px)>1, hxbk = - sum(px.*log2(fastchangem(px,1,0))); else hxbk = 0; end;
if length(py)>1, hybg = - sum(py.*log2(fastchangem(py,1,0))); else hybg = 0; end;

if (length(px)*length(py)==1)|(hxbk==0)|(hybg==0),
 p=0;
else
 p = sum(sum(pdf.*log2(fastchangem(pdf,1,0)./(py*px)))) / sqrt(hxbk * hybg);
end;

