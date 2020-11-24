function total_data_batch=SC_remove_batches( total_data,idx_samples, idx_batches)
% SC_remove_batches ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% Removes batch effects 

% INPUT
% total_data: expression matrix
% idx_samples: indexes of the conditions  (like wt1, wt2, ko1 ....)
% idx_samples: indexes of the batches 

% OUTPUT
% total_data_batch: expression matrix batch removed


% setting random seed
rng(1)

sum_ex = sum(total_data);
total_data_norm=total_data ./ repmat(sum_ex,length(total_data),1) * mean(sum_ex);

% CORRECTION=1 only in rare cases with unbalanced batch design
CORRECTION=0;

for i=1:length(idx_samples)
    for j=1:length(idx_batches)
        present(i,j)=any(intersect(idx_samples{i},idx_batches{j}));
    end
end

if CORRECTION
    good_pools=find(sum(present)==max(sum(present)))
    bad_conditions=find(sum(present')==max(sum(present')))
end

clear present i j

num_genes=length(total_data_norm(:,1));

% improvable maybe adding more percentiles or adjusting them to the
% effective size of the groups.
percentiles=[0:100];

% fixing batches in each condition at a time
for conditions=1:length(idx_samples)

    conditions
    
    all_indices=[];
    indexes={};
    for k=1:length(idx_batches)
        indexes{k}=intersect(idx_samples{conditions},idx_batches{k});
        all_indices = [all_indices indexes{k}'];
    end
    
    % removing empty pools, very important !
    indexes=indexes(cellfun(@any,indexes));
    
    
    for h=1:num_genes
        
        if nnz(total_data_norm(h,all_indices))
        
            vq = zeros(length(indexes),length(percentiles));
            tot = zeros(length(indexes),1);
            
            for k=1:length(indexes) % test pool by pool
                vq(k,:) = prctile(total_data_norm(h,indexes{k}),percentiles);
                tot(k,1) = length(indexes{k});
                [ dummy ix] = sort(total_data_norm(h,indexes{k}));
                nums=unique(dummy);
                for n=1:length(nums)
                    pos=find(dummy==nums(n));
                    if (length(pos)>1)    
                        values=ix(pos);
                        ix(pos)=values(randperm(length(pos)));
                    end
                end
                order{k}=ix;
            end
            
            
            if (CORRECTION & any(intersect(conditions,bad_conditions)) )
                cleaned=sum(vq(good_pools,:).*repmat(tot(good_pools),1,length(vq(1,:))))/sum(tot(good_pools));
            else
                cleaned=sum(vq.*repmat(tot,1,length(vq(1,:))))/sum(tot);
            end
            
            for k=1:length(indexes)
                dummy=indexes{k};
                total_data_norm(h,dummy(order{k})) = interp1( [1:(tot(k)-1)/100:tot(k)],cleaned,[1:tot(k)]);
            end
        
        end
        
    end
    
end

bar(sum(total_data_norm))
total_data_batch=(total_data_norm/mean(sum_ex)).*repmat(sum_ex,length(total_data_norm),1);
total_data_batch=round(total_data_batch);


end

