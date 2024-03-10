function [T] = py_GenKI(X, g, idx, wkdir, isdebug)

if nargin < 5, isdebug = true; end
if nargin < 4, wkdir = []; end

oldpth = pwd();
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, 'external', 'py_GenKI');

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

%prgfoldername = 'py_GenKI';
%narginchk(3, 3);
%assert(size(X, 1) == length(g));

%T = [];
%[pyok, wrkpth, x] = run.pycommon(prgfoldername);
%if ~pyok
%    error('GenKI (requires Python) has not been installed or set up correctly.')
%end

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


%     pw1=fileparts(mfilename('fullpath'));
%     wrkpth=fullfile(pw1,'external','py_GenKI');
%     cd(wrkpth);
%
%     fw = gui.gui_waitbar([],[],'Checking Python environment...');
%     x=pyenv;
%     try
%         pkg.i_add_conda_python_path;
%     catch
%
%     end
%     cmdlinestr=sprintf('"%s" "%s%srequire.py"', ...
%             x.Executable,wrkpth,filesep);
%     disp(cmdlinestr)
%     [status,cmdout]=system(cmdlinestr,'-echo');
%     if status~=0
%         gui.gui_waitbar(fw,true);
%         cd(oldpth);
%         waitfor(errordlg(sprintf('%s',cmdout)));
%         error('Python GenKI has not been installed properly.');
%     end
%
%     if isvalid(fw)
%         gui.gui_waitbar(fw,[],'Checking Python environment is complete');
%     end

try
    tmpfilelist = {'X.mat', 'g.txt', 'c.txt', ...
        'idx.mat', 'output.csv', fullfile('GRNs', 'pcNet_example.npz')};
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
    if issparse(X), X=full(X); end
    save('X.mat', '-v7.3', 'X');
    save('idx.mat', '-v7.3', 'idx');
    writematrix(g, 'g.txt');
    % writematrix(ones(size(X, 2), 1), 'c.txt');
    % disp('Input X g c written.');
catch ME
    if isvalid(fw)
         gui.gui_waitbar(fw, true);
    end
    errordlg(ME.message,'');
    return;
end
if isvalid(fw)
    gui.gui_waitbar(fw, [], [], 'Checking Python environment is complete');
    pause(0.5);
    gui.gui_waitbar(fw, [], [], 'Running GenKI...');
end


if isvalid(fw)
    gui.gui_waitbar(fw, [], [], 'Building pcnet_Source network...');
end
    A1 = sc_pcnetpar(sce.X(:, sce.c_cell_type_tx == celltype1));
    A1 = A1 ./ max(abs(A1(:)));
    A = ten.e_filtadjc(A1, 0.75, false);
    save('pcnet_Source.mat', 'A', '-v7.3');
    % disp('pcnet_Source.mat saved.');
if isvalid(fw)
    gui.gui_waitbar(fw, [], [], 'pcnet_Source.mat saved.');
end

codefullpath = fullfile(codepth,'script.py');
cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
disp(cmdlinestr)
[status] = system(cmdlinestr, '-echo');

if status == 0 && isvalid(fw)
    gui.gui_waitbar(fw, [], 'DataMapPlot is complete');
end



% fw=gui.gui_waitbar([],[],'Running GenKI...');

%    cmdlinestr=sprintf('"%s" "%s%sscript.py"', ...
%        x.Executable,wrkpth,filesep);
%    disp(cmdlinestr)
%    [status]=system(cmdlinestr,'-echo');

%[status] = run.pycommon2(x, wrkpth, prgfoldername);

% if isvalid(fw)
%     gui.gui_waitbar(fw,[],'Running GenKI is complete');
% end

if status == 0 && exist('output.csv', 'file')
    T = readtable('output.csv');
    T.Properties.VariableNames{1} = 'gene';
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
