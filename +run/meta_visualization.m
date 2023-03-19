
pw1=fileparts(mfilename('fullpath'));
if ~(ismcc || isdeployed)    
    pth=fullfile(pw1,'thirdparty','PHATE'); % for calling randmds.m
    addpath(pth);
    pth1=fullfile(pw1,'thirdparty','umapFileExchange');
    pth3=fullfile(pw1,'thirdparty','umapFileExchange','umap.jar');
    addpath(pth1);
    javaaddpath(pth3);
end

