function [ScoreV,T]=sc_cellcyclescoring(X,genelist)
% https://rdrr.io/cran/Seurat/src/R/utilities.R
% https://satijalab.org/seurat/v2.4/cell_cycle_vignette.html
% https://github.com/theislab/scanpy/blob/5db49c2612622b84b6379922204e0ba453d80ca2/scanpy/tools/_score_genes.py
% 
% # Trying here to match the Seurat approach in scoring cells.
% # Basically we need to compare genes against random genes in a matched
% # interval of expression.

% [ScoreV]=run.SeuratCellCycle(X,genelist);

X=sc_norm(X);
X=zscore(X);
X(X>10)=10;   % https://github.com/satijalab/seurat/issues/1166

exprv=sparse_nanmean(X,2);

nbins=25;
rng default
%[~,~,bin]=histcounts(log(20+exprv),nbins);
[bin]=discretize(log(10+exprv),nbins);

[~,sgenes,g2mgenes]=pkg.i_get_cellcyclegenes;
g=upper(genelist);
[~,idx1]=intersect(g,sgenes);
[~,idx2]=intersect(g,g2mgenes);

i1=i_randsmplbin(idx1,bin);
i2=i_randsmplbin(idx2,bin);

%i1=randsample(length(g),length(idx1));
%i2=randsample(length(g),length(idx2));

%%

score_s=sparse_nanmean(X(idx1,:))-sparse_nanmean(X(i1,:));
score_g=sparse_nanmean(X(idx2,:))-sparse_nanmean(X(i2,:));

ScoreV=string(repmat('G1',size(X,2),1));
score_S=score_s';
score_G2M=score_g';
C=[score_S,score_G2M];
[~,I]=max(C,[],2);
Cx=C>0;
i=sum(Cx,2)>0;
ScoreV(i&I==1)="S";
ScoreV(i&I==2)="G2M";
if nargout>1
    T=table(score_S,score_G2M,ScoreV);
end
% (score_s>score_g)&(sum(Cx,2)==2);  %'S'
% (score_s<score_g)&(sum(Cx,2)==2);  %'G2M'
% (Cx(:,1)>0)&(Cx(:,2)==0);          %'S'
% (Cx(:,1)==0)&(Cx(:,2)>0);          %'G2M'
end


function y=sparse_nanmean(X,dim)
    if nargin<2, dim=1; end
    if dim~=1, X=X'; end
    Z=isnan(X)|(X==0);
    z=sum(~Z);
    y=sum(X,'omitnan');
    y=y./z;
    if dim~=1, y=y'; end
end

function i=i_randsmplbin(idx,bin)
    i=[];
    binpool=bin(idx);
    ubinpool=unique(binpool);
    for k=1:length(ubinpool)
        tmpi=binpool==ubinpool(k);
        i=[i; randsample(find(bin==ubinpool(k)),sum(tmpi))];
    end
end
