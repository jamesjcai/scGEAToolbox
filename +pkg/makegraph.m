function G=makegraph(A,genelist)
    if isa(A,'digraph')||isa(A,'graph')
        G=A;
    else
        if issymmetric(A)
            G=graph(A,genelist,'omitselfloops');    
        else
            G=digraph(A,genelist,'omitselfloops');            
        end
    end
end
