function [y,z_new,b_new,cut_indices,I] = modifySpacingAllComps(y,z,b,x,mass,I,r,cut_indices,num_add,eo)
% This function makes sure that the spacing between the points on the
% curve(s) is appropriate. Re-spacing is done on all curves simultaneosly, 
% and at every function call, either only even or only odd numbered indices
% are considered, as specified by input parameter 'eo'. Additional points 
% may be added, as indicated by the input parameter 'num_add'.
%
% INPUTS: y-current points for curve
% mass-the masses of the data points
% I-the assignment of the data points to y (corresponding to shortest distance)
% x-the data points
% r-the threshold for modifying the spacing of points on y
% cut_indices-the last y indices for each component
% num_add-the net gain in number of points on the curve to be desired
% eo-whether to check the odd or even number indices: 0-even, 1-odd
% p- alpha_i^(p-1) * y_mass(i) should be constant
%
% OUTPUTS: y_new - the new points for the curve
% I_new - the estimated assignments for y_new (only the assignments for
% the two neighbors of added/removed points are re-evaluated)


assert(r>=1, 'The value of r should be at least 1.');
assert(eo==0 || eo==1, 'The value of eo should be either 0 or 1.');
[m,~] = size(y);

% compute mass at each point in y, and remove points with zero mass.
[y,z_new,b_new,I,y_mass,cut_indices] = remZeroMassPoints(y,z,b,mass,I,cut_indices);
[m1,d] = size(y);
k1 = m-m1;
m = m1;

% compute alpha(i): the relative length of the interval around y(i)
% extending to the midpoints of it's 2 neighbors.
w = zeros(m+2,d);
w(1,:) = y(1,:);
w(2:m+1,:) = y(:,:);
w(m+2,:) = y(m,:);

v = w(2:m+2,:)-w(1:m+1,:);
normv = sqrt(sum(v.^2,2));
normv(cut_indices+1) = 0;
total_length = sum(normv);

alpha = .5*(total_length/m)*(normv(2:m+1) + normv(1:m));
%compute alpha(i)^(p-1)*y_mass(i): want it to be roughly constant for each i
am = bsxfun(@times,alpha,y_mass);

am_mean = mean(am);

comp_lengs = [cut_indices;m]-[0;cut_indices];
comp_start_indices = [1;cut_indices+1];
nonsingleton_indices = setdiff(1:m,comp_start_indices(comp_lengs == 1));

start_pts = [1;cut_indices+1]; %the start and end indices for all components of y
end_pts = [cut_indices;m];     %to be checked in the loop below.
nonendpts = setdiff(1:m,[start_pts;end_pts]);

[am_sorted,sort_ind] = sort(am);
removals_bool = (am_sorted(1:floor(m/2))/am_mean < 1/r) & mod(sort_ind(1:floor(m/2)),2)==eo ...
    & ismember(sort_ind(1:floor(m/2)),nonendpts);
removals = sort_ind([removals_bool; boolean(zeros(ceil(m/2),1))]);

k = length(removals);
if ~isempty(num_add)
    k2 = min(m-k,k1+k+num_add);
    k2 = min(k2,floor(m/2));
else
    k2 = k;
end

split = sort_ind(m-k2+1:m);
split = sort(intersect(split,nonsingleton_indices),'descend');
k2 = length(split);

num_false_add = 0; %number of points which don't get added

for i=1:k2 
    
    if sum(split(i) == start_pts) > 0 %add after split(i)
        I_rel = I == split(i) | I==split(i)+1;
        t0 = .5;
        I_new = [];
        while sum(I_new==3) == 0 && t0>.05
            new_point = (1-t0)*y(split(i),:)+t0*y(split(i)+1,:);
            I_new = dsearchn([y(split(i),:);y(split(i)+1,:);new_point],x(I_rel,:));
            t0 = t0/2;
        end
        if sum(I_new==3) == 0
            t0 = .25;
            while sum(I_new==3) == 0 && t0>.05
                new_point = t0*y(split(i),:)+(1-t0)*y(split(i)+1,:);
                I_new = dsearchn([y(split(i),:);y(split(i)+1,:);new_point],x(I_rel,:));
                t0 = t0/2;
            end
        end
        if sum(I_new==3)>0
            add_ind = split(i)+1;
        else
            add_ind = [];
            num_false_add = num_false_add+1;
        end
    else %add before split(i)
        I_rel = I == split(i) | I==split(i)-1;
        t0 = .5;
        I_new = [];
        while sum(I_new==3) == 0 && t0>.05
            new_point = (1-t0)*y(split(i),:)+t0*y(split(i)-1,:);
            I_new = dsearchn([y(split(i),:);y(split(i)-1,:);new_point],x(I_rel,:));
            t0 = t0/2;
        end
        if sum(I_new==3) == 0
            t0 = .25;
            while sum(I_new==3) == 0 && t0>.05
                new_point = t0*y(split(i),:)+(1-t0)*y(split(i)-1,:);
                I_new = dsearchn([y(split(i),:);y(split(i)-1,:);new_point],x(I_rel,:));
                t0 = t0/2;
            end
        end
        if sum(I_new==3) > 0
            add_ind = split(i);
        else
            add_ind = [];
            num_false_add = num_false_add+1;
        end
    end
    start_pts(start_pts>=split(i)) = []; %to optimize checking condition
    end_pts(end_pts>=split(i)) = [];
    
    if ~isempty(add_ind)
        [y,z_new,b_new] = addPoints(y,z_new,b_new,new_point,add_ind,1);
        %need for re-indexing after changes
        removals = removals + double(removals>=add_ind);
        cut_indices = cut_indices + double(cut_indices>=add_ind);
        m = m+1;
    end
    
end

y_mass = zeros(m,1);
I = dsearchn(y,x)';
for i=1:m
    y_mass(i) = sum(mass(I==i));
end
nomass = find(y_mass == 0);
removals = sort([removals(1:max(length(removals)-num_false_add-length(nomass),0));nomass],'descend');
removals = sort(removals(1:end-num_false_add),'descend');
[y,z_new,b_new] = removePoints(y,z_new,b_new,removals);
k = length(removals);
for i=1:k
    cut_indices = cut_indices - double(cut_indices>=removals(i));
    I = I - double(I>=removals(i));
end

cut_indices(cut_indices<1)=[];

end

