function [T, Tf, OUT] = sc_tenifoldnet(X0, X1, genelist)
%scTenifoldNet pipeline for constructing scGRNs from scRNAseq data X0 and
%X1 and comparing two constructed networks
%
% [T,Tf,OUT]=sc_tenifoldnet(X0,X1,genelist)
%
% X0 and X1 are two gene-by-cell matrices
%
% use ten.sctenifoldnet for more input parameters.

import ten.*

if exist('sctenifoldnet', 'file') ~= 2
    error('Requires sctenifoldnet.m');
end

[T, A0, A1] = ten.sctenifoldnet(X0, X1, genelist);
glist = T.genelist(T.pAdjusted < 0.05);
run.web_Enrichr(glist);

Tf = e_fgsearun(T);
OUT = e_fgseanet(Tf, "JaccardCutoff", 0.5, "PlotNetwork", false, "ShowNotepad", false);

g = OUT{1, 1};
[A0sml, A1sml, glst] = gui.i_extractnet2compare(A0, A1, genelist, g);
sc_grnview2(A0sml, A1sml, glst);

end