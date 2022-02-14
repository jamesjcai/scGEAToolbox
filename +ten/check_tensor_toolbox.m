function check_tensor_toolbox

if exist(['@tensor',filesep,'tensor.m'],'file')~=2
    if ~(ismcc || isdeployed)
    try
        folder = fileparts(mfilename('fullpath'));
        a=strfind(folder,filesep);
        folder=extractBefore(folder,a(end)+1);
        pth=fullfile(folder,'+run','thirdparty','tensor_toolbox');
        addpath(pth);
    catch
            error(sprintf('TENSOR TOOLBOX is not installed.\nTo download, visit https://www.tensortoolbox.org/'));
    end
    end
end
if exist(['@tensor',filesep,'tensor.m'],'file')~=2
    disp('Visit https://www.tensortoolbox.org, to download TENSOR TOOLBOX.');
    error(sprintf('TENSOR TOOLBOX is not installed.\nTo download, visit https://www.tensortoolbox.org/'));
end
