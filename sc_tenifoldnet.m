function [T,Tf,OUT]=sc_tenifoldnet(X0,X1,genelist)

[T,A0,A1]=sctenifoldnet_m(X0,X1,genelist);
glist=T.genelist(T.pAdjusted<0.05);
run.Enrichr(glist);

Tf=e_fgsearun(T);
OUT=e_fgseanet(Tf,"JaccardCutoff",0.5,"PlotNetwork",false,"ShowNotepad",false);

g=OUT{1,1};
[A0sml,A1sml,glst]=gui.i_extractnetworks2compare(A0,A1,genelist,g);
sc_grnview2(A0sml,A1sml,glst);

end
