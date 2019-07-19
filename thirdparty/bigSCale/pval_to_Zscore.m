function [ z_score ] = pval_to_Zscore(pval,capovolgi)

x=[-100:0.005:100];

distribuzione=normcdf(x,0,1);

if (capovolgi)
    pval(pval>=0.5)=0.499999999;
end

[m n]=size(pval);


for i=1:m
for k=1:n
    
    if (pval(i,k)==1)
        z_score(i,k)=0; % modificato a luglio 2015
    else
        %z_score(i)=x(max(find(distribuzione<=pval(i))));
        z_score(i,k)=x(sum(distribuzione<=pval(i,k)));
    end
end

end


z_score(abs(z_score)<=0.005)=0;

%if (~iscell(pval))
% else
%     
%     
%     
% dummy_distribuzione=repmat(distribuzione,length(pval{1,1}),1);
%     
% for i=1:length(pval(:,1))
%     for j=1:length(pval(1,:))
%         j
%         pval_in_use=pval{i,j};
%         if (capovolgi)
%         pval_in_use(pval_in_use>=0.5)=0.499999999;
%         end
%         
%         test=dummy_distribuzione<=repmat(pval_in_use,1,length(dummy_distribuzione));
%         monitor=sum(test');
%         z_score{i,j}=x(monitor);
%         
%     end
% end
%     
% end


end

