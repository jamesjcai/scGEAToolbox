function [sout]=harmonypy(s,batchid,usepylib)

isdebug=false;

if nargin<3, usepylib=false; end
if nargin<2, error('[s]=run.harmonypy(s,batchid)'); end
oldpth=pwd();
pw1=fileparts(mfilename('fullpath'));
wrkpth=fullfile(pw1,'external','py_harmonypy');
cd(wrkpth);

tmpfilelist={'output.csv','input1.csv','input2.csv'};

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end


%if exist('./output.csv','file'), delete('./output.csv'); end
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
    
    try
        pkg.i_add_conda_python_path;
    catch        
    end
    cmdlinestr=sprintf('"%s" "%s%srequire.py"', ...
            x.Executable,wrkpth,filesep);
    disp(cmdlinestr)
    [status,cmdout]=system(cmdlinestr,'-echo');
    if status~=0
        cd(oldpth);
        sprintf('%s\n',cmdout);
        % waitfor(errordlg(sprintf('%s',cmdout)));
        error('harmony-py has not been installed properly.');
    end

    cmdlinestr=sprintf('"%s" "%s%sold_script.py"',x.Executable,wrkpth,filesep);
    disp(cmdlinestr)
    [status]=system(cmdlinestr);
    if status==0 && exist('output.csv','file')
        sout=readmatrix('output.csv');
    else
        sout=[];
    end
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end

