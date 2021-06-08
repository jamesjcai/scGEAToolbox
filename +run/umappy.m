function [sout]=umappy(X,usepylib)
if nargin<2, usepylib=false; end
if nargin<1, error('[s]=run.umappy(X,true)'); end
oldpth=pwd();
pw1=fileparts(mfilename('fullpath'));
wrkpth=fullfile(pw1,'thirdparty','harmony');
cd(wrkpth);

if exist('output.csv','file'), delete('output.csv'); end
%writematrix(X,'input1.csv');
writetable(array2table(X),'input1.csv');

% pyenv('Version','d:\\Miniconda3\\envs\\harmonypy\\python.exe')

if usepylib
    pd = py.importlib.import_module('pandas');
    um = py.importlib.import_module('umap');
    data_mat = pd.read_csv("input.csv");
    embedding = um.UMAP().fit_transform(data_mat);
    sout=np2mat(embedding);
else    
    x=pyenv;
    cmdlinestr=sprintf('"%s" "%s%sscript.py"',x.Executable,wrkpth,filesep);
    disp(cmdlinestr)
    [status]=system(cmdlinestr);
    if status==0 && exist('output.csv','file')
        sout=readmatrix('output.csv');
    else
        sout=[];
    end
end
if exist('input1.csv','file'), delete('input1.csv'); end
if exist('input2.csv','file'), delete('input2.csv'); end
if exist('output.csv','file'), delete('output.csv'); end
cd(oldpth);
end

