function check_tensor_toolbox

if exist(['@tensor',filesep,'tensor.m'],'file')~=2    
    if ~(ismcc || isdeployed)
    try
        folder = fileparts(mfilename('fullpath'));
        %a=strfind(folder,filesep);
        %folder=extractBefore(folder,a(end)+1);
        % pth=fullfile(folder,'+run','thirdparty','tensor_toolbox');
        pth=fullfile(folder,'..','tensor_toolbox');
        addpath(pth);
        disp(pth)
    catch
        error('Path to TENSOR TOOLBOX cannot be added.');
            % error(sprintf('TENSOR TOOLBOX is not installed.\nTo download, visit https://www.tensortoolbox.org/'));
    end
    end
end
if exist(['@tensor',filesep,'tensor.m'],'file')~=2
    disp('Visit https://www.tensortoolbox.org, to download TENSOR TOOLBOX.');
    error(sprintf('TENSOR TOOLBOX is not installed.\nTo download, visit https://www.tensortoolbox.org/'));
end

end

function i_downloadtensortoolbox
%     pw1=cdgea;
%     %wrkpth=fullfile(pw1,'tensor_toolbox2');
%     url='https://gitlab.com/tensors/tensor_toolbox/-/archive/v3.5/tensor_toolbox-v3.5.zip';
%     %unzip(url,wrkpth);
%     aaa=fullfile(tempdir,'aaa.zip');
%     websave(aaa,url);
%     unzip(aaa,pw1);
    
    url='https://gitlab.com/tensors/tensor_toolbox/-/archive/v3.5/tensor_toolbox-v3.5.zip';
    try
        unzip(url,tempdir);
    catch ME
        error(ME.message);
    end
    srcpth=fullfile(tempdir,'tensor_toolbox-v3.5');
    if exist(srcpth,'dir')
        pw1=cdgea;
        wrkpth=fullfile(pw1,'tensor_toolbox');
        [status,message] = ...
            movefile(srcpth,wrkpth);
        if status~=1
            error(message);
        end
    else
        error('check_tensor_toolbox error.');
    end
end
