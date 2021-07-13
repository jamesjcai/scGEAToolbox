oldpath=pwd();
cdgea;
cd ../
tic
disp('Updating scGEAToolbox...')
unzip('https://github.com/jamesjcai/scGEAToolbox/archive/master.zip');
addpath('./scGEAToolbox-master');
toc
if exist('sc_scatter.m','file')
    disp('scGEAToolbox updated!')
end
cd(oldpath);
