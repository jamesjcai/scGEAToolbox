function [ScoreV, T] = sc_cellcyclescore(X, g)
    % Score cell cycle phases

    % Define path to cell cycle gene list
    pw1 = fileparts(mfilename('fullpath'));
    wrkpth = fullfile(pw1, 'resources', 'CellScores', 'cellcyclegenes.xlsx');
    
    % Read gene table from the file
    T = readtable(wrkpth, ...
        'ReadVariableNames', true, 'FileType', 'spreadsheet', ...
        'Sheet', 'Regev_cell_cycle_genes');
    
    % Extract S and G2M phase genes
    sgenes = string(T.S);
    sgenes = sgenes(strlength(sgenes) > 0);
    g2mgenes = string(T.G2M);
    g2mgenes = g2mgenes(strlength(g2mgenes) > 0);
    
    % Calculate scores for S and G2M phases
    score_S = sc_cellscore_admdl(X, g, sgenes);
    score_G2M = sc_cellscore_admdl(X, g, g2mgenes);
    
    % Assign cell cycle phase based on scores
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
    
    % Return a table if requested
    if nargout > 1
        T = table(score_S, score_G2M, ScoreV);
    end
end   
