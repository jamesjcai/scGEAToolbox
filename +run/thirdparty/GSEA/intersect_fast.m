function [common,indA,indB] = intersect_fast(uA,uB)
% Find intersect between two datasets

uA = uA(:); uB = uB(:);

% Find matching entries
[sortuAuB,indSortuAuB] = sort([uA;uB]);
indInterAB = sortuAuB(1:end-1)==sortuAuB(2:end);

% Force the correct shape when indIcInterAB is 0.
if isequal(indInterAB, 0)                   
    indInterAB = zeros(0,1);
end

% Find common set
common = sortuAuB(indInterAB);                   

% Find indices if needed.
if nargout > 1
    indA = 1:length(uA);
    indB = 1:length(uB);
    indInterAB = find(indInterAB);
    ndx = indSortuAuB([indInterAB;indInterAB+1]);
    lenA = length(uA);
    logIndab = ndx > lenA;          % Find indices by index values in ndx.
    indA = indA(ndx(~logIndab));
    indB = indB(ndx(logIndab)-lenA);
end