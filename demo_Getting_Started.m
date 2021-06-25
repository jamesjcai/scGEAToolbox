% Learn to use SC_SCATTER
echo off
clc
echo on
%% DEMO_GETTING_STARTED - Using SC_SCATTER
%             
% scGEAToolbox contains SC_SCATTER for exploratory scRNAseq data 
% analysis. For this demonstration you will need to view both the
% command window and figure windows.

%% (1/4) Change working directory

cdgea;

pause % Press any key to continue...

%% (2/4) Load example data sets
% X - gene-by-cell expression matrix
% g - gene list
% s - coordinates of cell embedding

load example_data/testXgs.mat X g s

pause  % Press any key to continue...

%% (3/4) Make a SingleCellExperiment object SCE

sce=SingleCellExperiment(X,g,s);

pause  % Press any key to continue...

%% (4/4) Run SC_SCATTER with SCE

sc_scatter(sce);

%% End of DEMO_GETTING_STARTED
echo off
 