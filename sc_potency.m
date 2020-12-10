function [r]=sc_potency(X,genelist,speciesid)
%estimate differentiation potency of cells
%without the need to assume prior biological knowledge 
%such as marker expression or timepoint.
% CCAT (Correlation of Connectome And Transcriptome)
%https://github.com/aet21/SCENT
%
%see also: sc_cytotrace

if nargin<3, speciesid=1; end

genelist=upper(genelist);
if speciesid==1
    ppinetfile='Z:\Cailab\CCC_utilities\STRING\stringdb_human.mat';
else
    ppinetfile='Z:\Cailab\CCC_utilities\STRING\stringdb_mouse.mat';
end
load(ppinetfile,'G');
G.Edges.Weight=double(G.Edges.Weight>0);
GNodes=upper(string(table2array(G.Nodes)));
Gdegree=G.degree;
X=log2(X+1.1);

[~,i,j]=intersect(genelist,GNodes);
d=Gdegree(j);
X=X(i,:);
r=corr(X,d);  % Correlation of Connectome And Transcriptome
