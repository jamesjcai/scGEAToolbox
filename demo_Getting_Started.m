echo off
clc
echo on
%% DEMO_GETTING_STARTED   Demonstration of SC_SCATTER Application.
%              James J. Cai.

% The scGEAToolbox contains SC_SCATTER Application for exploratory 
% single-cell RNA-seq (scRNA-seq) data analysis. For this demonstration
% you will need to view both the command window and figure windows.

%% === (1/5) Change current directory ===

cdgea;

pause  % Press any key to continue...

%% === (2/5) Load example data sets ===

load example_data/testXgs.mat X g s
% X - gene-by-cell expression matrix
% g - gene list
% s - coordinates of cell embedding

load example_data/testSce.mat sce
% sce - SingleCellExperiment object/class

pause  % Press any key to continue...

%% === (3/5) Run SC_SCATTER ===

sc_scatter(X,g,s);

pause  % Press any key to continue...

%% === (4/5) Run SC_SCATTER_SCE ===

sc_scatter_sce(sce);

pause  % Press any key to continue...

%% === (5/5) Example of multidimensional view ===

gui.sc_multiembeddings

%% End of DEMO_GETTING_STARTED
echo off
 