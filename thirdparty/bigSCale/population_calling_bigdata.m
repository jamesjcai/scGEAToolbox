function [driving_genes table]= population_calling_bigdata( indices, indptr, data, OPZIONE_coreg,min_driving,min_coreg)


reps=round(length(indptr-1)/20000)
if reps>50
    reps=50;
end

if min(indptr)==0 & min(indices)==0
    disp('Fixing indexes');
    indptr=indptr+1;
    indices=indices+1;
end

num_genes=max(indices);
num_samples=length(indptr)-1;

SIZE=10000;

for j=1:reps
    
    subset=unique(randi(num_samples,1,SIZE*2));
    subset=subset(randperm(length(subset)));
    subset=subset(1:SIZE);

    total_data=zeros(num_genes,SIZE,'single');

    for k=1:SIZE
        cell=subset(k);
        total_data(indices( indptr(cell) : indptr(cell+1) - 1 ) , k )=data( indptr(cell) : indptr(cell+1)-1 );
    end
    
    [driving_genes{j}]= population_calling_v2( total_data,0,min_driving,min_coreg);
    close all;
end


table=zeros(num_genes,reps);
for k=1:length(driving_genes)
    table(driving_genes{k},k)=1;
end
if reps>1
    driving_genes=find( (sum(table')/reps)>0.25);
else
    driving_genes=driving_genes{1};
end

end


 