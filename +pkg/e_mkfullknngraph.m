[A]=sc_knngraph(s,4);
G=graph(A);
b=conncomp(G);
for k=1:max(b)-1
    for l=k:max(b)
        s(b==k,:)
        s(b==l,:)
        (k,l)
    end
end