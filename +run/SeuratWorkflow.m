function [sce]=SeuratWorkflow(X,genelist)
%Run cell cycle analysis using R package Seurat
%Seurat implements the method proposed by Tirosh et al.39 to score cells based on the averaged normalized expression of known markers of G1/S and G2/M.
%https://science.sciencemag.org/content/352/6282/189

if isa(X,'SingleCellExperiment') && isnumeric(genelist)
    sce=X;
    ndim=genelist;
else
    if nargin<2, error("[c]=run.SeuratCellCycle(X,genelist)"); end
    sce=SingleCellExperiment(X,genelist);
    ndim=2;
end
oldpth=pwd();
[isok,msg]=commoncheck_R('R_SeuratWorkflow');
if ~isok, error(msg); end
if exist('input.txt','file'), delete('input.txt'); end
if exist('tsneoutput.csv','file'), delete('tsneoutput.csv'); end
if exist('umapoutput.csv','file'), delete('umapoutput.csv'); end
if exist('activeidentoutput.csv','file'), delete('activeidentoutput.csv'); end
sc_writefile('input.txt',sce.X,sce.g);
if ndim==3
    RunRcode('script3d.R');
else    
    RunRcode('script.R');
end
if exist('tsneoutput.csv','file')
    T=readtable('tsneoutput.csv','ReadVariableNames',true);
    s=table2array(T(:,2:end));
    sce.struct_cell_embeddings.tsne=s;
    sce.s=s;
end
if exist('umapoutput.csv','file')
    T=readtable('umapoutput.csv','ReadVariableNames',true);
    s=table2array(T(:,2:end));
    sce.struct_cell_embeddings.umap=s;
end
if exist('activeidentoutput.csv','file')
    T=readtable('activeidentoutput.csv','ReadVariableNames',true);
    c=table2array(T(:,2:end));
    sce.c_cluster_id=c;
    sce.struct_cell_clusterings.seurat=c;
    sce.c=c;
end
if exist('input.txt','file'), delete('input.txt'); end
if exist('tsneoutput.csv','file'), delete('tsneoutput.csv'); end
if exist('umapoutput.csv','file'), delete('umapoutput.csv'); end
if exist('activeidentoutput.csv','file'), delete('activeidentoutput.csv'); end
cd(oldpth);
end