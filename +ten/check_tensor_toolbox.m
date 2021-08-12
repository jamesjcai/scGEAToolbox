function check_tensor_toolbox
if exist(['@tensor',filesep,'tensor.m'],'file')~=2
    disp('Visit https://www.tensortoolbox.org, to download TENSOR TOOLBOX.');
    error(sprintf('TENSOR TOOLBOX is not installed.\nTo download, visit https://www.tensortoolbox.org/'));
end