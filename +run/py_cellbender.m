function [output_h5] = py_cellbender(input_h5, wkdir, isdebug)

% PMID: 37550580
output_h5 = [];
prgfoldername = 'py_cellbender';

if nargin<2 || isempty(wkdir)
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(prgfoldername, preftagname);
    if isempty(wkdir), return; end
end
if nargin < 3, isdebug = true; end



oldpth = pwd();
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, '..',  'external', prgfoldername);

if isempty(wkdir) || ~isfolder(wkdir)
    cd(codepth);
else
    disp('Using working directory provided.');
    cd(wkdir);
end

% fw = gui.gui_waitbar([], [], 'Checking Python environment...');

x = pyenv;
try
    pkg.i_add_conda_python_path;
catch

end

codepth = pkg.i_normalizepath(codepth);

codefullpath = fullfile(codepth,'require.py');
cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);

disp(cmdlinestr)
[status, cmdout] = system(cmdlinestr, '-echo');
if status ~= 0
    cd(oldpth);
    %if isvalid(fw)
    %    gui.gui_waitbar(fw, true);
    %end
    error(cmdout);
end
if nargin < 1, input_h5 = []; end
if ~exist(input_h5,"file")
    [filenm, pathname] = uigetfile( ...
        {'*.h5;*.hdf5', 'HDF5 Files (*.h5)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Pick 10x Genomics H5 file(s)','MultiSelect','off');
    if isequal(filenm, 0), return; end            
    input_h5 = fullfile(pathname, filenm);
end

tmpfilelist = {'input.txt', 'output.h5'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
% input_h5 = """"+input_h5+"""";

writelines(input_h5,"input.txt");
%save('input.mat', '-v7.3', 'input_h5');
disp('Input file written.');

%if isvalid(fw)
%    gui.gui_waitbar(fw, [], [], 'Checking Python environment is complete');
%    pause(0.5);
%    gui.gui_waitbar(fw, [], [], 'Running CellBender...');
%end
% fw = gui.gui_waitbar([],[],'Running Scrublet...');

if canUseGPU
    codefullpath = fullfile(codepth,'script_gpu.py');
else    
    codefullpath = fullfile(codepth,'script.py');
end
pkg.i_addwd2script(codefullpath, wkdir, 'python');

cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
disp(cmdlinestr)
[status] = system(cmdlinestr, '-echo');

if status == 0 && exist('output_filtered.h5', 'file')
    output_h5 = fullfile(pwd, 'output_filtered.h5');    
end

%if status == 0 && isvalid(fw)
%    gui.gui_waitbar(fw, [], 'CellBender is complete');
%end

% if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);

end
