% Learn to use SCGEATOOL - Single-cell Gene Expression Analysis Tool
echo off
clc
echo on
%% DEMO_GETTING_STARTED - Using SCGEATOOL
%             
% The scGEATool provides a flexible interface where you can interactively 
% explore scRNAseq data and perform analysis. For this demonstration you 
% will need to view both the command window and figure windows.

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

%% (4/4) Run SCGEATOOL with SCE

scgeatool(sce);

%% End of DEMO_GETTING_STARTED
echo off
 