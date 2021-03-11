function sc_grnview2(A1,A2,genelist)

if nargin<3, error('USAGE: sc_grnview2(A1,A2,g)'); end
G1=pkg.makegraph(A1,genelist);
G2=pkg.makegraph(A2,genelist);
gui.i_doublegraphs(G1,G2);

end
