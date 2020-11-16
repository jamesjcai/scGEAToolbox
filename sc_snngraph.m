function [A]=sc_snngraph(s)
% https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
% As in PhenoGraph, we first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). This step is performed using the FindNeighbors function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs).
%
% https://github.com/liurui39660/SNNDPC/tree/MatlabImplementation/src/m
%{
A=[   0   1   1   1   0
   1   0   1   1   0
   1   1   0   1   1
   1   1   1   0   1
   0   0   1   1   0];
%}
A=sc_knngraph(s,5);

B=zeros(size(A));
n=size(A,1);
for i=1:n-1
    for j=i+1:n
        if A(i,j)>0
            B(i,j)=sum(sum(A([i j],:))>1);
        end
    end
end
B=B+B';
A=B;
% ref: https://slideplayer.com/slide/3373523/ with a typo

% http://mlwiki.org/index.php/SNN_Clustering#SSN_Clustering_Algorithm