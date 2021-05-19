echo off
clc
echo on
%% DEMO_GETTING_STARTED   Demonstration of SC_SCATTER Application.
%              James J. Cai.

% The scGEAToolbox contains SC_SCATTER Application for exploratory 
% single-cell RNA-seq (scRNA-seq) data analysis. For this demonstration
% you will need to view both the command window and figure windows.

%% === (1/4) Change current directory ===

cdgea;

pause  % Press any key to continue...

%% === (2/4) Load example data sets ===

load example_data/testXgs.mat X g s
% X - gene-by-cell expression matrix
% g - gene list
% s - coordinates of cell embedding

pause  % Press any key to continue...

%% === (3/4) Make a SingleCellExperiment object SCE ===

sce=SingleCellExperiment(X,g,s);

pause  % Press any key to continue...

%% === (4/4) Run SC_SCATTER with SCE ===

sc_scatter(sce);


%% End of DEMO_GETTING_STARTED
echo off
 