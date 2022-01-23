tic
disp('Installing scTenifoldNet...')
unzip('https://github.com/cailab-tamu/scTenifoldNet/archive/master.zip');
addpath('./scTenifoldNet-master/MATLAB');
toc

tic
disp('Installing scGEAToolbox...')
unzip('https://github.com/jamesjcai/scGEAToolbox/archive/master.zip');
addpath('./scGEAToolbox-master');
toc
savepath;
if exist('cdgea.m','file')
    disp('scGEAToolbox installed!')
end
% webread('https://api.github.com/repos/jamesjcai/scgeatoolbox')

%rootSettings = matlab.internal.getSettingsRoot;
%addonFolder  = rootSettings.matlab.addons.InstallationFolder.ActiveValue;
%executableFolder = fullfile(addonFolder,'Apps',name_of_your_app);