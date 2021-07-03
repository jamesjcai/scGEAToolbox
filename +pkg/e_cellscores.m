function [score]=e_cellscores(X,genelist,type)
if nargin<3, type="T_Cell_Exhaustion"; end
switch type
    case "T_Cell_Exhaustion"
        tgsPos=["CD69","PDCD1","TGFB1","CTLA4","SPN","LAG3"];
        %tgsPos=["CD44","LY6C","KLRG1","CTLA","ICOS","LAG3"];
        tgsNeg="";
    case "T_Cell_Cytotoxicity"
        
    case "Macrophage_Polarization_Index"
        tgsPos=["TGFB1"];
        tgsNeg=["Retnla"];
    otherwise
        error('Undefined')
end
[score]=sc_cellscore(X,genelist,tgsPos,tgsNeg);
end