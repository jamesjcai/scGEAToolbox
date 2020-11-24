function [ FC ] = SC_DE_fc( total_data, indici )

somma_ex=sum(total_data);
num_genes=length(total_data(:,1));

total_data_norm=total_data./repmat(somma_ex,num_genes,1)*mean(somma_ex);

f1=sum(total_data_norm(:,indici{1})'>0)./sum(total_data_norm(:,indici{1})'>=0);
f2=sum(total_data_norm(:,indici{2})'>0)./sum(total_data_norm(:,indici{2})'>=0);
    
e1=zeros(size(f1));
e2=zeros(size(f1));

for k=1:num_genes
	e1(k)=mean(nonzeros(total_data_norm(k,indici{1})));
	e2(k)=mean(nonzeros(total_data_norm(k,indici{2})));
        
	if isnan(e1(k))
        e1(k)=0;
    end
	if isnan(e2(k))
        e2(k)=0;
    end
        
	if f1(k)>f2(k)
        e2(k) = (f2(k)*e2(k)) / (f1(k)) ;
    else
        e1(k) = (f1(k)*e1(k)) / (f2(k)) ;
    end
    
end


FC=log2(e2'./e1');


end

