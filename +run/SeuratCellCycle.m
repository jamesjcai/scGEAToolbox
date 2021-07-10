function [c,T]=SeuratCellCycle(X,genelist)
%Run cell cycle analysis using R package Seurat
%Seurat implements the method proposed by Tirosh et al.39 to score cells based on the averaged normalized expression of known markers of G1/S and G2/M.
%https://science.sciencemag.org/content/352/6282/189
if nargin<2, error("[c]=run.SeuratCellCycle(X,genelist)"); end
oldpth=pwd();
[isok,msg]=commoncheck_R('R_SeuratCellCycle');
if ~isok, error(msg); end

if exist('output.csv','file')
    delete('output.csv');
end
sc_writefile('input.txt',X,upper(genelist));
pkg.RunRcode('script.R');


if exist('output.csv','file')
    T=readtable('output.csv','ReadVariableNames',true);
    c=string(T.Phase);
else
    c=[]; T=[];
end
if exist('input.txt','file'), delete('input.txt'); end
if exist('output.csv','file'), delete('output.csv'); end
cd(oldpth);
end