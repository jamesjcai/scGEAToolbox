function [nmad] = sc_qcmetrics(X)
    % SC_QCMETRICS Calculate quality control metrics for single-cell data
    %
    %   NMAD = SC_QCMETRICS(X)
    %
    %   Input:
    %       X - Gene expression matrix (genes x cells)
    %
    %   Output:
    %       NMAD - Number of median absolute deviations below the median
    %              of quality metrics
    %
    %   Description:
    %   This function calculates quality control metrics for single-cell data.
    %   The NMAD value represents the number of median absolute deviations
    %   below the median of two quality metrics:
    %   1. Log of library size
    %   2. Log of number of expressed genes
    %   
    %   A cell is removed if it falls below one of these metrics by more than
    %   the specified NMAD value. Higher NMAD values remove fewer cells,
    %   while lower values remove more cells.
    %
    %   Reference:
    %   https://alona.panglaodb.se/add.html

    % Calculate number of cells
    nc = size(X, 2);
    
    % Calculate percentage of cells expressing each gene
    pct = sum(X > 0, 2) ./ nc; % 1%
    
    % Calculate reads per cell
    reads_per_cell = sum(X, 2)'; % reads per cell
    
    % Calculate NMAD
    nmad = mad(log(sum(X, 2)) - log(sum(X > 0, 2)), 1);
end