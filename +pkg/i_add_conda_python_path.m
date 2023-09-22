x = pyenv;
py_root_useFromMATLAB = fileparts(x.Executable);
ENV = getenv('PATH');
ENV = strsplit(ENV, ';');
items_to_add_to_path = {; ...
    fullfile(py_root_useFromMATLAB, 'Library', 'mingw-w64', 'bin'); ...
    fullfile(py_root_useFromMATLAB, 'Library', 'usr', 'bin'); ...
    fullfile(py_root_useFromMATLAB, 'Library', 'bin'); ...
    fullfile(py_root_useFromMATLAB, 'Scripts'); ...
    };
ENV = [string(items_to_add_to_path(:)); string(ENV(:))];
ENV = unique(ENV, 'stable');
ENV = strjoin(ENV, ';');
setenv('PATH', ENV);

% module_to_load = 'gseapy';
% python_module_to_use = py.importlib.import_module(module_to_load);
% py.importlib.reload(python_module_to_use);

% https://www.mathworks.com/matlabcentral/answers/443558-matlab-crashes-when-using-conda-environment-other-than-base
