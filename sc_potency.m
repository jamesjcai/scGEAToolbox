function [r]=sc_potency(X,genelist)
%estimate differentiation potency of cells
%without the need to assume prior biological knowledge 
%such as marker expression or timepoint.
% CCAT (Correlation of Connectome And Transcriptome)
%https://github.com/aet21/SCENT

genelist=upper(genelist);
ppinetfile='Z:\Cailab\CCC_utilities\STRING\stringdb.mat';
load(ppinetfile,'G');
[~,i,j]=intersect(G.Nodes,genelist,"stable");

X=X(j,:);
genelist=genelist(j);
G=subgraph(G,i);
X=log2(X+1.1);
[~,idx]=ismember(genelist,G.Nodes);
% genelist=genelist(idx);
X=X(idx,:);
r=corr(X,G.degree);  % Correlation of Connectome And Transcriptome



