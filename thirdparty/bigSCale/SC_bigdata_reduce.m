function [ total_data , indices, indptr, data] = SC_bigdata_reduce( indices, indptr, data , genes )

%

indices = uint16(indices);
data = uint16(data);
genes = uint16(genes);

if min(indptr)==0 & min(indices)==0
    disp('Fixing indexes');
    indptr=indptr+1;
    indices=indices+1;
end

num_samples=length(indptr)-1



result=ismember(indices,genes);


detected=zeros(1,num_samples);

for k=1:length(indptr)-1
    detected(k)=sum(result(indptr(k) : indptr(k+1)-1));
end

if min(detected==0)
    error('Error unknown');
end

indices=indices(result);
[~, ~, indices]=unique(indices);
data=data(result);
indptr=cumsum([1 detected]);

total_data=zeros(numel(genes),num_samples ,'uint16');

%total_data=sparse(numel(genes),num_samples);

for k=1:num_samples
    total_data(indices( indptr(k) : indptr(k+1)-1 ),k)=data( indptr(k) : indptr(k+1)-1 );
end