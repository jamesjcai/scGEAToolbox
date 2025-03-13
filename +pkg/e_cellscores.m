function [score, T, posg] = e_cellscores(X, genelist, ...
    scoretypeid, methodid, showwaitbar)
% Calcute predefined cell scores (marker list in cellscores.xlsx)
%
% see also: SC_CELLSCORE_UCELL, SC_CELLSCORE_ADMDL, SC_CELLCYCLESCORING

% scoretypeid

if nargin < 5, showwaitbar = true; end
if nargin < 4, methodid = []; end
if nargin < 3, scoretypeid = 0; end
if nargin < 2, genelist = []; end
if nargin < 1, X = []; end

pw1 = fileparts(mfilename('fullpath'));

try
    cellscoresfile = fullfile(pw1, '..', 'resources', 'CellScores', ...
        'cellscores.xlsx');
    T = readtable(cellscoresfile, 'Sheet', 'Sheet1', ...
        'ReadVariableNames', true);
catch ME
    disp(ME.message);
    cellscoresfile = fullfile(pw1, '..', 'resources', 'CellScores', 'cellscores.txt');
    T = readtable(cellscoresfile, 'Delimiter', '\t', ...
        'ReadVariableNames', true);
end

%T=sortrows(T,"ScoreType");

if ischar(scoretypeid) || isstring(scoretypeid)
    idx = find(matches(T.ScoreType, scoretypeid, 'IgnoreCase', true));
elseif isnumeric(scoretypeid)
    idx = scoretypeid;
end
if isempty(idx)
    score = [];
    return;
end
idx = idx(1);    % in case there is duplicate in the score name list
if ~(idx <= height(T) && idx > 0 && idx == floor(idx))
    score = [];
    return;
end

scoretype = string(T.ScoreType(idx));

tgsPos = unique(strsplit(string(T.PositiveMarkers(idx)), ','));
tgsNeg = unique(strsplit(string(T.NegativeMarkers(idx)), ','));
if ~isempty(tgsNeg)
    tgsNeg = setdiff(tgsNeg, tgsPos);
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


%[score]=sc_cellscore_admdl(X,genelist,tgsPos,tgsNeg);
%[score]=sc_cellscore_ucell(X,genelist,tgsPos);

% if methodid==1
%     answer=gui.i_pickscoremethod(1); %'UCell [PMID:34285779]';
% elseif methodid==2
%     answer=gui.i_pickscoremethod(2); %'AddModuleScore/Seurat';
% else
%     answer = gui.i_pickscoremethod;
% end

answer = gui.i_pickscoremethod(methodid);
switch answer
    case 'AddModuleScore/Seurat'
        if showwaitbar, fw = gui.gui_waitbar; end
        try
            [score] = sc_cellscore_admdl(X, genelist, tgsPos, tgsNeg);
        catch ME
            if showwaitbar, gui.gui_waitbar(fw, true); end
            errordlg(ME.message);
            return;
        end
        if showwaitbar, gui.gui_waitbar(fw); end
    case 'UCell [PMID:34285779]'
        %[cs]=run.UCell(sce.X,sce.g,posg);

        if showwaitbar, fw = gui.gui_waitbar([], [], scoretype); end
        try
            [score] = sc_cellscore_ucell(X, genelist, tgsPos);
        catch ME
            if showwaitbar, gui.gui_waitbar(fw, true); end
            errordlg(ME.message);
            return;
        end
        if showwaitbar, gui.gui_waitbar(fw); end
    otherwise
        return;
end


% --------------------

posg = sort(tgsPos);

isexpressed = ismember(upper(posg), upper(genelist));
fprintf('\n=============\n%s (%s)\n-------------\n', 'Genes', scoretype);
for k = 1:length(posg)
    if isexpressed(k)
        fprintf('%s*, ', posg(k));
    else
        fprintf('%s, ', posg(k));
    end
    if mod(k, 10) == 0 || k == length(posg)
        fprintf('\n');
    end
end
fprintf('=============\n*Expressed genes (n = %d)\n', ...
    sum(isexpressed));

% fprintf('\n=============\n%s\n-------------\n','Marker Genes');
% for k=1:length(posg)
%     fprintf('%s\n',posg(k));
% end
% fprintf('=============\n');

end
