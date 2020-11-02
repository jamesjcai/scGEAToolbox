function [ ynew, znew, bnew ] = addPoints(y,z,b,points,locations,opt)
%ADDPOINTS adds new points to y, and adjusts z,b accordingly
%   if y = [2;4;23;65;26], points = [1;7], locations = [2,4], then the
%   output will be ynew = [2;1;4;7;23;65;26]
% z and b are also updated: entries for old values of z and b are kept, and
% the new entries are "reset".
% Use opt = 1 if points are being added to existing curve, opt = 0 if they
% are just added to singleton.

[m,d] = size(y);
[k,d1] = size(points);
assert(d==d1);

ynew = zeros(m+k,d);
ynew(locations,:) = points;
ynew(setdiff(1:m+k,locations),:) = y;

%initialize z,b
znew = zeros(m+k-1,d);
bnew = zeros(m+k-1,d);
znew(setdiff(1:m+k-1,locations),:) = [z;zeros(sum(locations==m+k),d)];
bnew(setdiff(1:m+k-1,locations),:) = [b;zeros(sum(locations==m+k),d)];

if opt == 0
    locs = union(locations,locations-1);
    locs = union(locs, locations+1);
    locs = locs(locs>0 & locs<m+k);
    % reset new entries
    znew(locs,:) = ynew(locs+1,:) - ynew(locs,:);
    bnew(locs,:) = 0;
else
    %need to modify z and b accordingly
    locs1 = locations-1;
    locs1 = locs1(locs1>0 & locs1<m+k-1);
    if ~isempty(locs1)
        d1 = sqrt(sum((ynew(locs1,:)-ynew(locs1+1,:)).^2,2))./sqrt(sum((ynew(locs1,:)-ynew(locs1+2,:)).^2,2));
        assert(d1<=1);
        znew(locs1,:) = bsxfun(@times,d1,znew(locs1,:));
        bnew(locs1,:) = bsxfun(@times,d1,bnew(locs1,:));
        locs2 = locations(locations>1 & locations<m+k);
        znew(locs2,:) = bsxfun(@times,(1-d1)/d1,znew(locs2-1,:));
        bnew(locs2,:) = bsxfun(@times,(1-d1)/d1,bnew(locs2-1,:));
    end
    locs3 = locations(locations == 1);
    if ~isempty(locs3)
        znew(locs3,:) = ynew(locs3+1,:) - ynew(locs3,:);
        bnew(locs3,:) = zeros(length(locs3),d);
    end
    locs4 = locations(locations == m+k);
    if ~isempty(locs4)
        znew(locs4-1,:) = ynew(locs4,:) - ynew(locs4-1,:);
        bnew(locs4-1,:) = zeros(length(locs4),d);
    end
end

end