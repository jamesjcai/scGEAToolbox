function [ ynew,znew,bnew ] = removePoints( y,z,b,locations )
%REMOVEPOINTS removes points from y.
%   This function deletes the points in y in the specified indices, and z
%   and b are updated appropriately. The 'locations' is a column vector
%   that specifies which indices in y to remove.

[m,~] = size(y);
ynew = y;
znew = z;
bnew = b;
if ~isempty(locations)
    adjust_locs = locations(locations>1 & locations<m);
    znew(adjust_locs,:) = z(adjust_locs,:)+z(adjust_locs-1,:);
    bnew(adjust_locs,:) = b(adjust_locs,:)+b(adjust_locs-1,:);
    
    rem_locs = adjust_locs-1;
    if sum(locations == 1)>0
        rem_locs = [rem_locs; min(setdiff(1:m-1,rem_locs))];
    end
    if sum(locations== m)>0
        rem_locs = [rem_locs; max(setdiff(1:m-1,rem_locs))];
    end
    rem_locs = rem_locs(rem_locs>0);
    znew(rem_locs,:) = [];
    bnew(rem_locs,:) = [];
    
    ynew(locations,:) = [];
end

end

