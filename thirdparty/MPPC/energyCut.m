function [ y_new,z_new,b_new,I_new,cut_indices_new ] = energyCut(x,mass,y,z,b,I,cut_indices,lambda,lambda2)
% This function cuts (removes) edges based on the energy.
% Single edges longer than lambda2 are removed which are not energetically
% advantageous, leaving the adjacent points (possibly as singletons). If an
% edge is shorter than lambda2, this function may remove an entire sequence
% of edges with total length > lambda2. In this case all vertices in the 
% sequence are removed as well.
% Note: 'I' may be passed as empty, in which case it will be computed inside
% of 'computeProjs'.

[m,~] = size(y);
k = length(cut_indices);
I_partial = computeProjs(x,y,I,cut_indices);
I_new = round(I_partial);
if ~isempty(I) && ~isequal(I_new,I)
    warning('I is not equal to I_new in energyCut');
end

v = y(2:m,:)-y(1:m-1,:);
normv = sqrt(sum(v.^2,2)); % edge lengths
normv(cut_indices) = 0;

cut_indices_new = union(cut_indices,find(normv>4*lambda2)); %cut edges longer than 4*lambda2
rem_ind=[];
i=0;
j=1;
while i<m-1
    i = i+1;
    if j<=k && i==cut_indices(j)
        j = j+1;
        continue;
    end
    curr = i;
    ind = I_partial>i & I_partial<i+1; %indices of data that project onto edge i
    I_partial0 = mod(I_partial(ind),1);
    edelta_fid = ((min(I_partial0,ones(size(I_partial0))-I_partial0)*normv(i)).^2)*mass(ind)';
    if lambda*normv(i) - lambda*lambda2 > edelta_fid % remove single edge
        cut_indices_new = union(cut_indices_new,i-length(rem_ind));
    else %check sequence of edges
        if normv(i) < lambda2 
            tot_len = normv(i);
            while curr<m-1
                curr = curr+1;
                if j<=k && curr==cut_indices(j)
                    j = j+1;
                    i = curr;
                    break;
                end
                 tot_len = tot_len + normv(curr);
                if tot_len>lambda2
                    ind0 = I_partial>i & I_partial<curr+1; %indices of data which project onto path considered for removal
                    tot_fid = mass(ind0)*sum((getYVals(y,I_partial(ind0))-x(ind0,:)).^2,2);
                    rest_ind = [i,curr+1]; % data will project here in most cases
                    [~,projdist] = dsearchn(y(rest_ind,:),x(ind0,:));
                    cut_fid = mass(ind0)*(projdist.^2); %corresponding fidelity if cut happens  
                    if tot_fid + lambda*tot_len > cut_fid + lambda*lambda2 %then remove entire path
                        rest_ind = setdiff(1:m,[rem_ind',i+1:curr]); % compute exact projections for new I
                        [I0,~] = dsearchn(y(rest_ind,:),x(ind0,:));
                        I_partial(ind0) = rest_ind(I0);
                        I_new = I_new - (curr-i)*(I_new>curr-length(rem_ind));
                        I_new(ind0) = I0'; 
                        new_cut_ind = i-length(rem_ind);
                        cut_indices_new = union(cut_indices_new,new_cut_ind); %cut entire sequence
                        rem_ind = union(rem_ind,i+1:curr); %singletons along this path will be removed
                        if ~iscolumn(rem_ind)
                            rem_ind = rem_ind';
                        end
                        cut_indices_new = cut_indices_new - length(i+1:curr)*(cut_indices_new>new_cut_ind);
                        i = curr; %advance index to next after cuts
                        break;
                    else
                        tot_len = tot_len - normv(i) - normv(curr);
                        i = i+1;
                        curr = curr - 1; %since it gets incremented later. want it to stay constant. 
                    end
                end
            end
        end
    end
end
cut_indices_new = cut_indices_new(cut_indices_new<m-length(rem_ind));
s = size(cut_indices_new);
if max(s)>1 && s(1)<s(2)
    cut_indices_new = cut_indices_new';
end
if isempty(cut_indices_new)
    cut_indices_new = []; %to avoid Matlab's concatenation warnings
end
[y_new,z_new,b_new] = removePoints(y,z,b,rem_ind);
end


function yvals = getYVals(y,partial_ind)
[~,d] = size(y);
partials = mod(partial_ind,1)';
yvals = repmat((1-partials),1,d).*y(floor(partial_ind),:) + repmat(partials,1,d).*y(ceil(partial_ind),:);

end
