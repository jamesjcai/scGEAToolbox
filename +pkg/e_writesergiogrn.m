close all
rng(244)
A=rand(5);
A=A-diag(diag(A));
B=A>0.55;
g=digraph(B);
plot(g)

d=indegree(g);
A=full(adjacency(g));
a1=triu(A);
a2=tril(A);


for k=1:size(g.Nodes,1)
    s='';
    %[d(k) sum(A(:,k))]
    f=find(A(:,k));
    f2=2*ones(size(f));
    f3=0.5*ones(size(f));
    s=sprintf('%d, %d, %s%s%s',k,d(k), ...
        sprintf('%d, ',f), ...
        sprintf('%f, ',f3), ...
        sprintf('%d, ',f2));
    disp(strip(s));

end
