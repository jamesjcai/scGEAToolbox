function [sout]=harmonypy(s,batchid,usepylib)
if nargin<3, usepylib=false; end
if nargin<2, error('[s]=run.harmonypy(s,batchid)'); end
oldpth=pwd();
pw1=fileparts(mfilename('fullpath'));
wrkpth=fullfile(pw1,'thirdparty','harmony');
cd(wrkpth);

if exist('output.csv','file'), delete('output.csv'); end
%writematrix(s,'input1.csv');
writetable(array2table(s),'input1.csv');
batchidx=matlab.lang.makeValidName(string(batchid));
writetable(table(batchidx),'input2.csv','QuoteStrings',true);
% pyenv('Version','d:\\Miniconda3\\envs\\harmonypy\\python.exe')

if usepylib    
    pd = py.importlib.import_module('pandas');
    np = py.importlib.import_module('numpy');
    hm = py.importlib.import_module('harmonypy');
    data_mat=pd.read_csv("input1.csv");
    data_mat=np.array(data_mat);
    meta_data = pd.read_csv("input2.csv");
    vars_use = py.list({py.str('batchidx')});
    ho = hm.run_harmony(data_mat, meta_data, vars_use);
    sout=np2mat(ho.Z_corr.T);
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

