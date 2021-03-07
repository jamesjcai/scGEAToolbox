function sc_grnview2(A1,A2,g)

if nargin<3, error('USAGE: sc_grnview2(A1,A2,g)'); end
G1=makeg(A1,g);
G2=makeg(A2,g);
gui.i_doublegraphs(G1,G2);

    function G=makeg(B,g)
        if isa(B,'digraph')||isa(B,'graph')
            G=B;
        else
            if issymmetric(B)
                G=graph(B,g,'omitselfloops');    
            else
                G=digraph(B,g,'omitselfloops');
            end
        end
    end
end
