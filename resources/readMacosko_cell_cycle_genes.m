T=readtable('Macosko_cell_cycle_genes.txt');
mm=(string(T.IG1_S));
mm(all(strcmp(mm,""),2),:) = [];
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4481139/
% 