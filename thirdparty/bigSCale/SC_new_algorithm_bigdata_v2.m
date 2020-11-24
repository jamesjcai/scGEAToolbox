function [N N_pct difference] = SC_new_algorithm_bigdata_v2( indices, indptr, data , edges)

SIZE=2000;
N_pct=[];
counts=1;
margin=10;
stop=95; % when xx % has an error within margin stop
max_cycles=5;

if min(indptr)==0 & min(indices)==0
    disp('Fixing indexes');
    indptr=indptr+1;
    indices=indices+1;
end

num_genes=max(indices)
num_samples=length(indptr)-1


while counts<max_cycles
    
    subset=unique(randi(num_samples,1,SIZE*2));
    subset=subset(randperm(length(subset)));
    subset=subset(1:SIZE);

    total_data=zeros(num_genes,SIZE,'double');

    for k=1:SIZE
        cell=subset(k);
        total_data(indices( indptr(cell) : indptr(cell+1) - 1 ) , k )=data( indptr(cell) : indptr(cell+1)-1 );
    end
   
    
    [N(:,:,counts) N_pct(:,:,counts) ] = SC_new_algorithm( total_data , edges, 0 );

    if counts>1
        
        [ r c z]=size(N_pct);
        
        for i=1:r
            for j=1:c
                media(i,j)=mean(N_pct(i,j,:));
                s_dev(i,j)=std(N_pct(i,j,:));
            end
        end
        
        sem=s_dev/sqrt(counts);
        difference=sem./media*100;
        
        sprintf('Cycle %g, found %g/%g (%g %%) values with the margin of %g %%',counts,nnz(difference<margin),numel(difference),nnz(difference<margin)/numel(difference)*100,margin)
        pause(5);
        
    end

    

    if counts>1 & nnz(difference<margin)/numel(difference)*100 > stop
        break;    
    end
    
    counts=counts+1
    
end

disp('Computing final values');
N=sum(N,3);
clear N_pct
% calculating percentiles
for k=1:length(N)
    for j=1:length(N)
        N_pct(k,j) = sum(N(k,j:end))/sum(N(k,:));
    end
end
N_pct(isnan(N_pct))=1;

end
