function [ScoreV, T] = sc_cellcyclescore(X, g)
    % Score cell cycle phases
    
    % [~,sgenes,g2mgenes]=pkg.i_get_cellcyclegenes;
    pw1 = fileparts(mfilename('fullpath'));
    wrkpth = fullfile(pw1, 'resources', 'cellcyclegenes.xlsx');
    
    T = readtable(wrkpth, ...
        'ReadVariableNames', true, 'FileType', 'spreadsheet', ...
        'Sheet', 'Regev_cell_cycle_genes');
    sgenes = string(T.S);
    sgenes = sgenes(strlength(sgenes) > 0);
    g2mgenes = string(T.G2M);
    g2mgenes = g2mgenes(strlength(g2mgenes) > 0);
    
    score_S = sc_cellscore_admdl(X, g, sgenes);
    score_G2M = sc_cellscore_admdl(X, g, g2mgenes);
    
    if all(isnan(score_S)) || all(isnan(score_G2M)) 
        ScoreV = string(repmat('unknown', size(X, 2), 1));
    else
        ScoreV = string(repmat('G1', size(X, 2), 1));
        C = [score_S, score_G2M];
        [~, I] = max(C, [], 2);
        Cx = C > 0;
        i = sum(Cx, 2) > 0;
        ScoreV(i & I == 1) = "S";
        ScoreV(i & I == 2) = "G2M";
    end
    if nargout > 1
        T = table(score_S, score_G2M, ScoreV);
    end
end
    
    
    % function [g,sgenes,g2mgenes]=i_get_cellcyclegenes
    % %Get cell-cycle genes
    %
    % % https://satijalab.org/seurat/v2.4/cell_cycle_vignette.html
    % % Also read in a list of cell cycle markers, from Tirosh et al, 2015
    % %cc.genes <- readLines(con = "~/Downloads/seurat_resources/regev_lab_cell_cycle_genes.txt")
    %
    % folder=fileparts(mfilename('fullpath'));
    % a=strfind(folder,filesep);
    % folder=extractBefore(folder,a(end)+1);
    % wrkpth=fullfile(folder,'resources',filesep);
    % fname=[wrkpth 'regev_lab_cell_cycle_genes.txt'];
    % g=string(readcell(fname));
    %
    % sgenes=g(1:43);
    % g2mgenes=g(44:97);
    %
    % return;
    % %https://science.sciencemag.org/content/352/6282/189
    %
    % % https://www.frontiersin.org/articles/10.3389/fgene.2017.00001/full
    %
    % %{
    % s1="CDK2, MKI67, RB1, E2F1, HIST1H2AE, CCNB1, CCNB1, CBX3, ND1, MKI67, GAPDH, CCNB1, CCNB2";
    % g1=unique(strsplit(s1,", "));
    % s2="ATAD2, BORA, BTG3, CASP8AP2, CCDC134, CCNA2, CCNB1, CCNB2, CCNE1, CCNE2, CDC25A, CDC25B, CDC25C, CDC5L, CDC7, CDCA3, CDCA4, CDCA7, CDK1, CDK19, CDK4, CDKN2C, CDKN2D, CDKN3, CKS1B, CKS2, CRLF3, DBF4B, DLGAP5, E2F1, E2F7, E2F8, FAM83D, FBXO5, FOXM1, FZR1, GMNN, LIN54, LIN9, MASTL, MELK, MYBL2, NPAT, ODF2, PA2G4, PBK, PIMREG, PKMYT1, PRR11, RBL1, STIL, TCF19, TFDP1, TICRR, TOE1, TRIM28, TTF2, UBE2C, UBE2S, UHRF2, USP37, WEE1, CDK2, AURKA, KIF14, DTL, GTSE1, PLK1";
    % g2=unique(strsplit(s2,", "));
    % g=unique([g1 g2]);
    %
    % return;
    %
    % options = weboptions('Timeout',21);
    % fname=tempname;
    % websave(fname,'https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv',options);
    % % t=readtable('a.txt','PreserveVariableNames',true);
    % warning off
    % t=readtable(fname,'Range','A:B');
    % warning on
    % g=string(t.ApprovedSymbol);
    % delete(fname);
    % %}
    %
    % % https://academic.oup.com/jmcb/article/11/8/703/5188008