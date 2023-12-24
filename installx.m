% tic
% disp('Installing scTenifoldNet...')
% unzip('https://github.com/cailab-tamu/scTenifoldNet/archive/master.zip');
% addpath('./scTenifoldNet-master/MATLAB');
% toc
tic
disp('Installing scGEAToolbox...')
unzip('https://github.com/jamesjcai/scGEAToolbox/archive/main.zip');
addpath('./scGEAToolbox-main');
toc
if exist('scgeatool.m','file')
    disp('scGEAToolbox installed!')
end
savepath(fullfile(userpath,'pathdef.m'));
% webread('https://api.github.com/repos/jamesjcai/scgeatoolbox')
