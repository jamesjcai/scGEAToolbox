function [score]=e_cellscores(X,genelist,type)
if nargin<3, type="T_Cell_Exhaustion"; end
switch type
    case "T_Cell_Exhaustion"
        tgsPos=["CD44","LY6C","KLRG1","CTLA","ICOS","LAG3"];
        tgsNeg=["IL2","TNF"];
    otherwise
        error('Undefined')
end
[score]=sc_cellscore(X,genelist,tgsPos,tgsNeg);