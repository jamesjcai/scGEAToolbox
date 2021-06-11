addpath("gaimc\graphs");
addpath("gaimc");
run("gspbox\gsp_start.m");
run("unlocbox\init_unlocbox.m");
%{
cd gspbox
gsp_start;
cd ../unlockbox
init_unlocbox;
cd ..
%}