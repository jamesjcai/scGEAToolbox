addpath("Z:\Cailab\Imputation_Test\gaimc\graphs");
addpath("Z:\Cailab\Imputation_Test\gaimc");
run("Z:\Cailab\Imputation_Test\gspbox\gsp_start.m");
run("Z:\Cailab\Imputation_Test\unlocbox\init_unlocbox.m");
%{
cd gspbox
gsp_start;
cd ../unlockbox
init_unlocbox;
cd ..
%}