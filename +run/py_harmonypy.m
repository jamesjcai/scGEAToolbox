function [sout] = py_harmonypy(s, batchid, wkdir, isdebug)
arguments
    s(:, :) {mustBeNumeric}
    batchid(1, :) {mustBePositive, mustBeInteger}
    wkdir = []
    isdebug = true
end
%if nargin < 4, isdebug = true; end
%if nargin < 3, wkdir = []; end

% prgfoldername = 'py_harmonypy';

oldpth = pwd();
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, '..', 'external', 'py_harmonypy');

if isempty(wkdir) || ~isfolder(wkdir)
    cd(codepth);
else
    disp('Using working directory provided.');
    cd(wkdir);
end

fw = gui.gui_waitbar([], [], 'Checking Python environment...');

x = pyenv;
try
    pkg.i_add_conda_python_path;
catch

end

codefullpath = fullfile(codepth,'require.py');
%cmdlinestr = sprintf('"%s" "%s%srequire.py"', ...
%    x.Executable, codepth, filesep);
cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);

disp(cmdlinestr)
[status, cmdout] = system(cmdlinestr, '-echo');
if status ~= 0
    cd(oldpth);
    if isvalid(fw)
        gui.gui_waitbar(fw, true);
    end
    error(cmdout);
end

sout = [];
% isdebug = false;
% 
% %if nargin<3, usepylib=false; end
% %if nargin<2, error('[s]=run.harmonypy(s,batchid)'); end
% 
% oldpth = pwd();
% [pyok, wrkpth, x] = run.pycommon(prgfoldername);
% if ~pyok, return; end

tmpfilelist = {'input.mat', 'output.mat'};
%  if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

%pw1=fileparts(mfilename('fullpath'));
%wrkpth=fullfile(pw1,'external','py_harmonypy');
%cd(wrkpth);
%tmpfilelist={'output.csv','input1.csv','input2.csv'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

if issparse(s), s = full(s); end
save('input.mat', '-v7.3', 's', 'batchid');
disp('Input file written.');



if isvalid(fw)
    gui.gui_waitbar(fw, [], [], 'Checking Python environment is complete');
    %pause(0.5);
    gui.gui_waitbar(fw, [], [], 'Running Harmonypy...');
end
codefullpath = fullfile(codepth,'script.py');

pkg.i_addwd2script(codefullpath, wkdir, 'python');
%     fw=gui.gui_waitbar([],[],'Running harmonypy...');
%     cmdlinestr=sprintf('"%s" "%s%sscript.py"', ...
%         x.Executable,wrkpth,filesep);
%     disp(cmdlinestr)
%     [status]=system(cmdlinestr,'-echo');
%     if isvalid(fw)
%         gui.gui_waitbar(fw,[],'Running harmonypy is complete');
%     end

cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
disp(cmdlinestr)
[status] = system(cmdlinestr, '-echo');


if status == 0 && exist('output.mat', 'file')
    load("output.mat", "sout")
end

if status == 0 && isvalid(fw)
    gui.gui_waitbar(fw, [], 'Harmonypy is complete');
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);

%if exist('./output.csv','file'), delete('./output.csv'); end
%writematrix(s,'input1.csv');
% writetable(array2table(s),'input1.csv');
% batchidx=matlab.lang.makeValidName(string(batchid));
% writetable(table(batchidx),'input2.csv','QuoteStrings',true);
% % pyenv('Version','d:\\Miniconda3\\envs\\harmonypy\\python.exe')
%
% if usepylib
%     pd = py.importlib.import_module('pandas');
%     np = py.importlib.import_module('numpy');
%     hm = py.importlib.import_module('harmonypy');
%     data_mat=pd.read_csv("input1.csv");
%     data_mat=np.array(data_mat);
%     meta_data = pd.read_csv("input2.csv");
%     vars_use = py.list({py.str('batchidx')});
%     ho = hm.run_harmony(data_mat, meta_data, vars_use);
%     sout=np2mat(ho.Z_corr.T);
% else
%     x=pyenv;
%
%     try
%         pkg.i_add_conda_python_path;
%     catch
%     end
%     cmdlinestr=sprintf('"%s" "%s%srequire.py"', ...
%             x.Executable,wrkpth,filesep);
%     disp(cmdlinestr)
%     [status,cmdout]=system(cmdlinestr,'-echo');
%     if status~=0
%         cd(oldpth);
%         sprintf('%s\n',cmdout);
%         % waitfor(errordlg(sprintf('%s',cmdout)));
%         error('harmony-py has not been installed properly.');
%     end
%
%     cmdlinestr=sprintf('"%s" "%s%sscript.py"',x.Executable,wrkpth,filesep);
%     disp(cmdlinestr)
%     [status]=system(cmdlinestr);
%     if status==0 && exist('output.csv','file')
%         sout=readmatrix('output.csv');
%     else
%         sout=[];
%     end
% end

end
