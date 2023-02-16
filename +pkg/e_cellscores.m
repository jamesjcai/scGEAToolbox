function [score,T,posg]=e_cellscores(X,genelist,typeid)
% Calcute predefined cell scores (marker list in cellscores.txt)
%
% see also: SC_CELLSCORE_UCELL, RUN.UCELL, SC_CELLSCORE, SC_CELLCYCLESCORING

if nargin<3, typeid=0; end
if nargin<2, genelist=[]; end
if nargin<1, X=[]; end

pw1=fileparts(mfilename('fullpath'));
cellscoresfile=fullfile(pw1,'..','resources','cellscores.txt');

T=readtable(cellscoresfile,'Delimiter','\t',...
    'ReadVariableNames',true);
%T=sortrows(T,"ScoreType");

if ischar(typeid) || isstring(typeid)
    idx=find(matches(T.ScoreType,typeid,'IgnoreCase',true));
elseif isnumeric(typeid)
    idx=typeid;
end
if isempty(idx)
    score=[];
    return;
end    
if ~(idx<=size(T,1) && idx>0 && idx == floor(idx))
    score=[];
    return;
end

tgsPos=unique(strsplit(string(T.PositiveMarkers(idx)),','));
tgsNeg=unique(strsplit(string(T.NegativeMarkers(idx)),','));
if ~isempty(tgsNeg)
    tgsNeg=setdiff(tgsNeg,tgsPos);
end

%{
if nargin<3, type="T_Cell_Exhaustion"; end
% https://carmonalab.github.io/UCell/UCell_vignette_TILstates.html#unsupervised-clustering
tgsNeg="";

switch type
    case "T_Cell_Exhaustion"
        tgsPos=["CD69","PDCD1","TGFB1","CTLA4","SPN","LAG3"];
        %tgsPos=["CD44","LY6C","KLRG1","CTLA","ICOS","LAG3"];
    case "T_Cell_Cytotoxicity"
    
    case "Macrophage_Polarization_Index"
        tgsPos=["TGFB1"];
        tgsNeg=["Retnla"];
    
	case "Fetal_epithelial_progenitor"
        tgsPos=["BEX3","STMN1","SOX4","LDHB","SKP1","SNRPE","ID3","SRP9","GSTP1","SRP14"];
    case "Macrophage"
        tgsPos=["CTSB","C1QB","LAPTM5","TYROBP","PSAP","C1QA","HLA-DRA","CTSD","NPC2","FCER1G"]; 
    case "B_cell_Plasmocyte"
        tgsPos=["JCHAIN","IGHA1","SSR4","MZB1","IGKC","IGHA2","HERPUD1","DERL3","SEC11C","FKBP11"];
    case "Fibroblast"
        tgsPos=["C1S","TIMP2","COL6A3","SEMA3C","MMP2","GSN","IGFBP6","MFAP4","COL6A1","PLAC9"];
    case "Fasciculata_cell"
        tgsPos=["PEBP1","STAR","RARRES2", "CLU","CYP21A2","CYP17A1","AKR1B1","NOV","TPD52L1", "EPHX1"];
    case "T_cell"
        tgsPos=["CD3D","CD3E","CD3G","CD4","CD2","CD7","TRAC","TRBC1","LAT"];        
    otherwise
        error('Undefined')
end
%}
[score]=sc_cellscore(X,genelist,tgsPos,tgsNeg);


posg=sort(tgsPos);

    isexpressed=ismember(upper(posg),upper(genelist));
    fprintf('\n=============\n%s\n-------------\n','Genes');
    for k=1:length(posg)
        if isexpressed(k)
           fprintf('%s\t*\n',posg(k));
        else
            fprintf('%s\t\n',posg(k));
        end
    end
    fprintf('=============\n*expressed genes (n=%d)\n', ...
        sum(isexpressed));

% fprintf('\n=============\n%s\n-------------\n','Marker Genes');
% for k=1:length(posg)
%     fprintf('%s\n',posg(k));
% end
% fprintf('=============\n');

end
