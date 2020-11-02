function [ num_add,avg_turn_deg ] = numPointsAdd( y, cut_indices, max_avg_turn, max_points, lambda2 )
%NUMPOINTSADD Determines the number of points to add to the existing curve.
%   The number of points to add will be > 0 as long as the average degree
%   of the turning angles is > the specified max_avg_turn deg and the current
%   number of points is less than max_points. At most 40%
%   of the current number of points on the curve will be added, in order to
%   prevent too many point being added simultaneously.
%   max_avg_turn is in degrees.

m = length(y(:,1));
avg_turn_deg = [];
if isempty(max_points) || m < max_points
    
    v = y(2:m,:)-y(1:m-1,:);
    normv = sqrt(sum(v.^2,2));
    num_add1 = 0;
    num_add2 = 0;
    if ~isempty(lambda2)
        normv(cut_indices) = 0;
        avg_len = sum(normv)/(m-1-length(cut_indices));
        if avg_len > lambda2/2
            num_add1 = max(1,floor(m * min(.4, 2*avg_len/lambda2))); 
        end
    end
    if ~isempty(max_avg_turn)
        vdots = sum(bsxfun(@times,v(2:m-1,:),v(1:m-2,:)),2);
        ta_deg = acosd(bsxfun(@rdivide,vdots,bsxfun(@times,normv(1:m-2),normv(2:m-1))));
        comp_ends = union(cut_indices(cut_indices>1)-1,cut_indices(cut_indices<m-1));
        ta_deg(comp_ends) = 0;
    
        avg_turn_deg = sum(ta_deg)/(m - 2 - length(comp_ends));
        
        if avg_turn_deg > max_avg_turn
            num_add2 = max(1,floor(m * min(.4, sqrt(avg_turn_deg/max_avg_turn - 1))));
        end
        start_ind = [1;cut_indices-1];
        end_ind = [cut_indices;m];
        num_add2 = min(num_add2+sum(end_ind == start_ind + 1),floor(.4*m)); %In case there are components with length 2 (no turning angle)
    end
    num_add = max(num_add1,num_add2);
    if ~isempty(max_points)
        num_add = min(max_points-m,num_add);
    end
else
    num_add = 0;
end

end

