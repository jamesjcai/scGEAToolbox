function [  somma_ex ] = SC_bigdata_somma_ex( indices, indptr, data )
%
somma_ex=zeros(1,numel(indptr)-1);

if min(indptr)==0 & min(indices)==0
    disp('Fixo gli indici');
    indptr=indptr+1;
    indices=indices+1;
end

for k=1:length(indptr)-1
    somma_ex(k)=sum(data(indptr(k) : indptr(k+1)-1));
end

% for k=1:length(indptr)-1
%     factor=somma_ex(k)/mean(somma_ex);
%     data(indptr(k) : indptr(k+1)-1)=data(indptr(k) : indptr(k+1)-1)/factor;
% end
