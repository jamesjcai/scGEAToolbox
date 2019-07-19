%

if min(indptr)==0 & min(indices)==0
    disp('Fixo gli indici');
    indptr=indptr+1;
    indices=indices+1;
end

num_genes=length(gene_names);
dummy=zeros(num_genes,1,'uint32');% !!! from 16!!!!!

indices_out=zeros(round(length(indices)),1,'uint16');
data_out=zeros(round(length(indices)),1,'uint32'); % !!! from 16!!!!!
indptr_out=zeros(length(blocks),1,'int64');

conta=1;

for k=1:length(blocks)
    k
	dummy(:)=0; 
    
	for j=1:nnz(blocks(k,:))
        cell=blocks(k,j);
        posizioni=indices(indptr(cell) : indptr(cell+1)-1);
        valori=data(indptr(cell) : indptr(cell+1)-1);
        dummy(posizioni)=dummy(posizioni)+valori;
    end
    
    idx=find(dummy);
    indices_out(conta:conta+length(idx)-1)=idx;
    data_out(conta:conta+length(idx)-1)=dummy(idx);
    indptr_out(k)=length(idx);
    conta=conta+length(idx);
end

indices_out=indices_out(1:conta-1);
data_out=data_out(1:conta-1);
indptr_out=cumsum([1 ; indptr_out]);