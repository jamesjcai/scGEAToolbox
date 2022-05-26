function sc_grnview2(A1,A2,nodenames,figname)
%GRN network viewer (two networks) 
    if nargin<4, figname=''; end
    if nargin<3, error('USAGE: sc_grnview2(A1,A2,g)'); end    
    G1=pkg.makegraph(A1,nodenames);
    G2=pkg.makegraph(A2,nodenames);
    gui.i_doublegraphs(G1,G2,figname);
end
