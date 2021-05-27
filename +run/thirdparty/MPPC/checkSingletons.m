function [y_new,z_new,b_new,cut_indices_new] = checkSingletons(x,mass,I,y,z,b,cut_indices,lambda,lambda2)
%checkSingletons adds and removes singletons to decrease MPPC energy

[m,d] = size(y);
cut_indices0 = cut_indices;
cut_indices_new = cut_indices;
start_indices = [1;cut_indices+1];
comp_lengs = [cut_indices;m] - start_indices + 1;
singletons = start_indices(comp_lengs==1);

k = length(singletons);
add_indices = zeros(k+m,1); %indices of the new singletons in y_new
add_points = zeros(k+m,d);
rem = zeros(k,1);
num_rem = 0;
num_add = 0;
hold on;
for i=1:k % loop through all existing singletons
    x_ind = ismember(I,singletons(i));
    n0 = sum(x_ind);
    fidelity = mass(x_ind)*sum((x(x_ind,:)-repmat(y(singletons(i),:),n0,1)).^2,2);
    y_rest_ind = setdiff(1:m,[singletons(i);rem]);
    if ~isempty(y_rest_ind)
        [I_xind,dists] = dsearchn(y(y_rest_ind,:),x(x_ind,:));
        fidelity_new = mass(x_ind)*(dists.^2);
    else
        fidelity_new = inf;
    end
    if fidelity_new - fidelity < lambda*lambda2 %then remove singleton
        rem(i) = singletons(i);
        I(x_ind) = y_rest_ind(I_xind);
        %scatter(y(rem(i),1),y(rem(i),2),'black','fill');
        num_rem = num_rem+1;
        ind = find(cut_indices==singletons(i));
        if isempty(ind) % singleton is last point in y
            cut_indices(end) = [];
            cut_indices_new(end) = [];
        else
            cut_indices(ind) = [];
            cut_indices_new = [cut_indices_new(1:ind-1);cut_indices_new(ind+1:end)-1];
        end
    else
        if lambda*lambda2 < fidelity % check whether to add singleton
            x1 = x(x_ind,:);
            x1 = x1(1,:);
            perturb = .5*(x1 - y(singletons(i),:));
            y0 = [y(singletons(i),:) + perturb; y(singletons(i),:) - perturb];
            [y0,~,dists] = EM(y0,x(x_ind,:),mass(x_ind),3);
            new_fid = mass(x_ind)*(dists.^2);
            if fidelity - new_fid > lambda*lambda2 %then add singleton
                num_add = num_add+1;
                new_ind = singletons(i) + num_add - num_rem;
                add_indices(i) = new_ind;
                if singletons(i) ~= m
                    ind = find(cut_indices==singletons(i));
                    cut_indices = [cut_indices(1:ind);0;cut_indices(ind+1:end)];
                    cut_indices_new = cut_indices_new + (cut_indices_new>=new_ind);
                    cut_indices_new = union(cut_indices_new,new_ind);
                else
                    cut_indices_new = union(cut_indices_new,new_ind-1);
                end
                y(singletons(i),:) = y0(1,:);
                add_points(i,:) = y0(2,:);
                %scatter(add_points(i,1),add_points(i,2),'red','fill');
            end
        else
            proj_dists = sqrt(sum((repmat(y(singletons(i),:),n0,1)-x(x_ind,:)).^2,2));
            tot_proj_mass = sum(mass(x_ind));
            av_proj_dist = mass(x_ind)*proj_dists/tot_proj_mass;
            
            diff = av_proj_dist - lambda/tot_proj_mass;
            if diff > 0  % if average proj distance > sqrt(lambda/est_density), add point so singleton can grow
                [max_dist,loc_max_ind] = max(proj_dists);
                max_ind = find(x_ind);
                max_ind = max_ind(loc_max_ind);
                perturb = .1*diff*(x(max_ind,:)-y(singletons(i),:))/max_dist;
                len = 2*norm(perturb);
                ydouble = [y(singletons(i),:)+perturb;y(singletons(i),:)-perturb];
                [~,new_dists] = dsearchn(ydouble,x(x_ind,:));
                new_fid = mass(x_ind)*(new_dists.^2);
                if new_fid + lambda*len < fidelity %check that energy indeed decreases
                    num_add = num_add+1;
                    new_ind = singletons(i) + num_add - num_rem;
                    add_indices(i) = new_ind;
                    cut_indices_new = cut_indices_new + (cut_indices_new>=new_ind-1);
                    y(singletons(i),:) = y(singletons(i),:) + perturb;
                    add_points(i,:) = y(singletons(i),:) - perturb;
                    %scatter(add_points(i,1),add_points(i,2),'red','fill');
                end
            end
        end
    end
end

% Now check non-singleton points along curve(s)
if isrow(cut_indices_new), cut_indices_new = cut_indices_new'; end
cut_indices_new = [cut_indices_new;zeros(m,1)];
num_add = 0;
proj_mass = zeros(m,1);
x_mean = zeros(m,d);
del_fidel = zeros(m,1);
del_leng = zeros(m,1);
j=1;
for i=1:m
    if isempty(cut_indices0) || i~=cut_indices0(j)
        proj_mass(i) = sum(mass(I==i));
        x_mean(i,:) = mass(I==i)*x(I==i,:)/proj_mass(i);
        del_fidel(i) = proj_mass(i)*sum((x_mean(i,:)-y(i,:)).^2,2);
        del_leng(i) = norm(y(i,:)-y(min(i+1,m),:)) + norm(y(i,:)-y(max(i-1,1),:)) - norm(y(min(i+1,m),:)-y(max(i-1,1),:));
        del_leng(singletons) = -inf;
        if lambda*lambda2 < del_fidel(i) + lambda*del_leng % then it is advantageous to add singleton
            add_points(m+i,:) = x_mean(i,:);
            add_indices(m+i) = m+1+num_add;
            cut_indices_new(m+i) = m+num_add;
            num_add = num_add+1;
            scatter(add_points(m+i,1),add_points(m+i,2),'red','fill');
        end
    else
        j = min(j+1,length(cut_indices0));
    end
end
hold off;

add_ind_bool = add_indices~=0;
add_indices = add_indices(add_ind_bool);
add_points = add_points(add_ind_bool,:);
cut_indices_new = cut_indices_new(cut_indices_new>0);

rem = rem(rem~=0);
if ~isempty(rem)
    [y,z,b] = removePoints(y,z,b,rem);
end
y_new = y;
z_new = z;
b_new = b;
if ~isempty(add_indices)
    [m,d] = size(y);
    k = length(add_indices);
    y_new = zeros(m+k,d); z_new = zeros(m+k-1,d); b_new = zeros(m+k-1,d);
    y_new(add_indices,:) = add_points;
    y_new(setdiff(1:m+k,add_indices),:) = y;
    add_indices1 = add_indices(add_indices>1);
    z_new(add_indices1-1,:) = zeros(length(add_indices1),d);
    b_new(add_indices1-1,:) = zeros(length(add_indices1),d);
    if length(add_indices) ~= length(add_indices1)
        z_new(1,:) = zeros(1,d);
        b_new(1,:) = zeros(1,d);
        z_new(setdiff(m+k-1,[1;add_indices-1]),:) = z;
        b_new(setdiff(m+k-1,[1;add_indices-1]),:) = b;
    else
        z_new(setdiff(1:m+k-1,add_indices-1),:) = z;
        b_new(setdiff(1:m+k-1,add_indices-1),:) = b;
    end
end

if isempty(cut_indices_new)
    cut_indices_new = []; %to avoid Matlab's concatenation warnings
end
end


function [y,I,dists] = EM(y0,x,mass,num_iter)
%   Helper EM function for checking when to add a singleton in the
%   neighborhood of another. y0 has two points, x is the data projecting to
%   them. num_iter is the number of desired iterations.

y = y0;
for i=1:num_iter
    I = dsearchn(y,x);
    I1 = I==1;
    y(1,:) = mass(I1)*x(I1,:)/sum(mass(I1));
    y(2,:) = mass(~I1)*x(~I1,:)/sum(mass(~I1));
end
[I,dists] = dsearchn(y,x);
I = I';
end

