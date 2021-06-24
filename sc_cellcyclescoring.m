function [ScoreV,T]=sc_cellcyclescoring(X,genelist)

[~,sgenes,g2mgenes]=pkg.i_get_cellcyclegenes;
score_S=sc_cellscore(X,genelist,sgenes);
score_G2M=sc_cellscore(X,genelist,g2mgenes);

ScoreV=string(repmat('G1',size(X,2),1));
C=[score_S,score_G2M];
[~,I]=max(C,[],2);
Cx=C>0;
i=sum(Cx,2)>0;
ScoreV(i&I==1)="S";
ScoreV(i&I==2)="G2M";
if nargout>1
    T=table(score_S,score_G2M,ScoreV);
end
end

