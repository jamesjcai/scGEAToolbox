function sc_grnview2(A1,A2,nodenamelist,figname)
    if nargin<4, figname=''; end
    if nargin<3, error('USAGE: sc_grnview2(A1,A2,g)'); end
    G1=pkg.makegraph(A1,nodenamelist);
    G2=pkg.makegraph(A2,nodenamelist);
    gui.i_doublegraphs(G1,G2,figname);
end
