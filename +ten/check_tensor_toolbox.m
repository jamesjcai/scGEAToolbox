function check_tensor_toolbox

if exist(['@tensor', filesep, 'tensor.m'], 'file') == 2
    return;
end

if ~(ismcc || isdeployed)
    if ispref('scgeatoolbox', 'tensor_toolbox_path')
        pth = getpref('scgeatoolbox', 'tensor_toolbox_path');
        if isfolder(pth)
            addpath(pth);
        end
    end
end

if exist(['@tensor', filesep, 'tensor.m'], 'file') ~= 2
    error('scgeatoolbox:tensortoolbox:notfound', ...
        'TENSOR TOOLBOX is not installed. Use the scTenifold menu to install it.');
end

end
