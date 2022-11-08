function [sout]=umappy(X,usepylib)
if nargin<2, usepylib=false; end
if nargin<1, error('[s]=run.umappy(X,true)'); end
oldpth=pwd();
pw1=fileparts(mfilename('fullpath'));
wrkpth=fullfile(pw1,'external','py_umappy');
cd(wrkpth);

if exist('./input.csv','file'), delete('./input.csv'); end
if exist('./output.csv','file'), delete('./output.csv'); end
X=sc_norm(X);
X=log(X+1);
if issparse(X), X=full(X); end
writematrix(X.','input.csv');
% writetable(array2table(X'),'input.csv','WriteVariableNames',false);

% pyenv('Version','d:\\Miniconda3\\envs\\harmonypy\\python.exe')

if usepylib
    pd = py.importlib.import_module('pandas');
    um = py.importlib.import_module('umap');
    data_mat = pd.read_csv("input.csv");
    embedding = um.UMAP().fit_transform(data_mat);
    sout=np2mat(embedding);
else
    x=pyenv;
    pkg.i_add_conda_python_path;
    cmdlinestr=sprintf('"%s" "%s%sscript.py"',x.Executable,wrkpth,filesep);
    disp(cmdlinestr)
    [status]=system(cmdlinestr);
    if status==0 && exist('output.csv','file')
        sout=readmatrix('output.csv');
    else
        sout=[];
    end
end
if exist('./input.csv','file'), delete('./input.csv'); end
if exist('./output.csv','file'), delete('./output.csv'); end
cd(oldpth);
end

