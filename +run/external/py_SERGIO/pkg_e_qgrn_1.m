load example_grn.mat
load output X
[M]=permn([0 1],size(X,1));  
y=0+(X>0).';
[~,idx]=ismember(y,M,'rows');
t=tabulate(idx);
n=zeros(size(M,1),1);
n(t(:,1))=t(:,2);

figure;
subplot(2,2,1)
bar(n)
for k=1:size(M,1)
    txt{k}=sprintf('%d',M(k,:));
end
txt=string(txt);
set(gca,'XTick',1:size(M,1));
set(gca,'XTickLabel',txt);
ylabel('# of cells');
xlabel('Expression pattern');


% ---------------------------------------------------------
f0=sum(y>0)./size(X,2);   % per gene initial activate freq.
pt=n./size(X,2);          % target cell state freq.    
% ---------------------------------------------------------


