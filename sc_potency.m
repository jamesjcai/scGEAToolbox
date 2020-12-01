function [r]=sc_potency(X,genelist)
%estimate differentiation potency of cells
%without the need to assume prior biological knowledge 
%such as marker expression or timepoint.
% CCAT (Correlation of Connectome And Transcriptome)
%https://github.com/aet21/SCENT

genelist=upper(genelist);
ppinetfile='Z:\Cailab\CCC_utilities\STRING\stringdb.mat';
load(ppinetfile,'G');
GNodes=string(table2array(G.Nodes));
Gdegree=G.degree;
X=log2(X+1.1);

[y,idx]=ismember(GNodes,genelist);
idx=idx(y);
% genelist=genelist(idx);
d=Gdegree(idx);
X=X(idx,:);

r=corr(X,d);  % Correlation of Connectome And Transcriptome



