function [T] = py_GSEApy_enr(genelist, backgroundlist, wkdir, ...
    showbarplot, isdebug)

if nargin < 5, isdebug = true; end
if nargin < 4, showbarplot = true; end
if nargin < 2, backgroundlist = []; end


extprogname = 'py_GSEApy_enr';
if nargin < 3 || isempty(wkdir)
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
    if isempty(wkdir), return; end
end

    oldpth = pwd();
    pw1 = fileparts(mfilename('fullpath'));
    codepth = fullfile(pw1, 'external', extprogname);
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


tmpfilelist = {'input.txt', 'background.txt',...
    'GO_Biological_Process_2023.Human.enrichr.reports.txt',...
    'GO_Molecular_Function_2023.Human.enrichr.reports.txt'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

writematrix(genelist, 'input.txt');
if ~isempty(backgroundlist)
    writematrix(backgroundlist, 'background.txt');
end
disp('Input files written.');

if isvalid(fw)
    gui.gui_waitbar(fw, [], [], 'Checking Python environment is complete');
    pause(0.5);
    gui.gui_waitbar(fw, [], [], 'Running GSEApy Enrichr...');
end
% fw = gui.gui_waitbar([],[],'Running DataMapPlot...');


    codefullpath = fullfile(codepth,'script.py');
    pkg.i_addwd2script(codefullpath, wkdir, 'python');
    cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
    disp(cmdlinestr)
    [status] = system(cmdlinestr, '-echo');


if status == 0 && isvalid(fw)

    gui.gui_waitbar(fw, [], 'GSEApy Enrichr is complete');
end

T1 = readtable('GO_Biological_Process_2023.Human.enrichr.reports.txt');
T2 = readtable('GO_Molecular_Function_2023.Human.enrichr.reports.txt');
T = [T1; T2];

if showbarplot
    figure; imshow(imread('GO_Biological_Process_2023.Human.enrichr.reports.png'))
    figure; imshow(imread('GO_Molecular_Function_2023.Human.enrichr.reports.png'))
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
