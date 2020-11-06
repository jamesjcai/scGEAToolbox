function [ y,cut_indices ] = connectCompsE(x,y,mass,I,cut_indices,lambda1,lambda2)
% This function connects curves based on the change in the continuous MPPC
% energy. The function uses a greedy approach to connecting curves 
% corresponding the largest reduction in energy. 

m = length(y(:,1));
k0 = length(cut_indices)+1;
k=k0;

start_ind = union(1,cut_indices+1);
if iscolumn(start_ind), start_ind = start_ind'; end
end_ind = union(cut_indices,m);
if iscolumn(end_ind), end_ind = end_ind'; end
leafs = [start_ind,end_ind];
assert(length(leafs)==2*k);


proj_indices = ismember(I,leafs);
proj = I(proj_indices); %projections of data that project to leafs
proj_x = x(proj_indices,:); %data that project to leafs
proj_mass = mass(proj_indices); %mass of data that project to leafs

D = computeDeltaE(y(leafs,:),leafs,proj,proj_x,proj_mass,lambda1,lambda2);
D(D==0) = inf;
D(sub2ind(size(D),1:k,k+1:2*k)) = inf;
D(sub2ind(size(D),k+1:2*k,1:k)) = inf;

D1 = reshape(D,4*k^2,1);
[vals,ind] = sort(D1);
[ind1,ind2] = ind2sub([2*k,2*k],ind);
num_add = min(sum(vals<0),k);

used_ends = [];
singletons = start_ind == end_ind;

for i0=1:num_add
    
    i = ind1(i0); j = ind2(i0);
    
    if isempty(cut_indices) || sum(ismember([i,j],used_ends))>0 || (...
            (ismember(mod(i-1,k)+1,mod(used_ends-1,k)+1)) && ...
            (ismember(mod(j-1,k)+1,mod(used_ends-1,k)+1)) && ...
            sum(singletons(mod([i,j]-1,k)+1))>0)
        continue;
    end
    used_ends = union(used_ends,[i,j]);
    
    if leafs(j)<leafs(i)
        temp = i;
        i = j;
        j = temp;
    end
    
    se = [[1;cut_indices+1],[cut_indices;m]];
    [i1,i2] = ind2sub(size(se),find(se==leafs(i)));
    [j1,j2] = ind2sub(size(se),find(se==leafs(j)));
    
    if i1==j1
        continue; % connecting would create loop
    end
    
    if length(i1)>1, i1 = i1(1); i2 = i2(1); end
    if length(j1)>1, j1 = j1(1); j2 = j2(1); end
    
    if i2==1  % connect to beginning of first component
        st_ind1 = se(i1,i2);
        e_ind1 = se(i1,2);
        comp1_ind = st_ind1:e_ind1;
        y1 = y(comp1_ind,:);
        y1 = flipud(y1);
        leng1 = 1+e_ind1-st_ind1;
        leafs(leafs==e_ind1) = st_ind1;
    else % connect to end of first component
        st_ind1 = se(i1,1);
        e_ind1 = se(i1,i2);
        comp1_ind = st_ind1:e_ind1;
        y1 = y(comp1_ind,:);
        leng1 = 1+e_ind1-st_ind1;
    end
    
    if j2==2  % connect to end of second component
        st_ind2 = se(j1,1);
        e_ind2 = se(j1,j2);
        comp2_ind = st_ind2:e_ind2;
        y2 = y(comp2_ind,:);
        y2 = flipud(y2);
        leng2 = 1+e_ind2-st_ind2;
        end_ind = leafs == st_ind2;
    else % connect to beginning of second component
        st_ind2 = se(j1,j2);
        e_ind2 = se(j1,2);
        comp2_ind = st_ind2:e_ind2;
        y2 = y(comp2_ind,:);
        leng2 = 1+e_ind2-st_ind2;
        end_ind = leafs == e_ind2;
    end
    
    y = [y(1:st_ind1-1,:);y1;y2;y(e_ind1+1:st_ind2-1,:);y(e_ind2+1:end,:)];
    if j1>length(cut_indices)
        j1 = j1-1;
    end
    cut_indices(j1) = [];
    cut_indices = cut_indices + leng2*(cut_indices>=leafs(i) & cut_indices<leafs(j));
    
    leafs = leafs+leng2*(leafs>e_ind1 & leafs < st_ind2);
    leafs(end_ind) = st_ind1+leng1+leng2-1;
end

end




function delta_E = computeDeltaE(ends,indices,proj,proj_x,proj_mass,lambda1,lambda2)
% Computes the change of energy by adding the edges (ends(i,:),end(k+1,:)).
% ends - consists of the endpoints which when connected to y0 form the
% edges of interest. Note that y0 is the last row of ends.
% indices - the indices corresponding to 'ends' in y.
% proj_x - the data that project to 'y0' or one of the positions in 'ends'.
% proj - the indices specifying to where the data in 'proj_x' project.
% proj(j) = i means that proj_x(j,:) projects to ends(i,:).
% proj_mass - the mass corresponding to the data 'proj_x'.

k = length(ends(:,1));
assert(mod(k,2) == 0)
k = k/2;
delta_E = zeros(2*k,2*k); % change in the energy

for i=1:2*k
    y0 = ends(i,:);
    ind1 = ismember(proj,indices(i));
    x1 = proj_x(ind1,:);
    n1 = length(x1(:,1));
    
    for j=[i+1:min(i+k-1,2*k),i+k+1:2*k]
        edge = [y0;ends(j,:)];
        edge_length = sqrt(sum((edge(1,:)-edge(2,:)).^2,2));
        if edge_length > 2*lambda2
            delta_E(i,j) = inf; % ensures we won't connect in this case
        else
            ind2 = ismember(proj,indices(j));
            x2 = proj_x(ind2,:);
            n2 = length(x2(:,2));
            edge_masses = [proj_mass(ind1),proj_mass(ind2)];
            % See how much data project to the interior of the edge
            I_partial = computeProjs([x1;x2],edge,[ones(1,n1),2*ones(1,n2)],[]);
            rel_ind = I_partial>1 & I_partial<2;
            I_partial = mod(I_partial(rel_ind),1);
            delta_fid = -((min(I_partial,ones(size(I_partial))-I_partial)*edge_length).^2)*edge_masses(rel_ind)';
            delta_E(i,j) = delta_fid + lambda1*edge_length - lambda1*lambda2;
        end
        delta_E(j,i) = inf;
    end
end

end

