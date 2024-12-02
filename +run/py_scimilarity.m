function [s] = py_scimilarity(sce, modeldir, wkdir, isdebug)
s = [];
if nargin < 4, isdebug = true; end
if nargin < 3, wkdir = tempdir; end
if nargin < 2
    modeldir = selectFolder;    
end
if isempty(modeldir) || ~exist(modeldir, 'dir')
    error('Model folder does not exist or is invalid.');
end

oldpth = pwd();
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, 'external', 'py_scimilarity');

if isempty(wkdir) || ~isfolder(wkdir)
    cd(codepth);
else
    disp('Using working directory provided.');
    cd(wkdir);
end
winopen(wkdir);

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
else 
    disp('Code requirement check is done.')
end

try
    pkg.i_deletefiles({'input.h5ad', 'output.h5ad'});
    tmpfilelist = {'Xnorm.mat', 'X.mat', 'g.csv', 'c.csv', ...
        'input.h5ad', 'output.h5ad'};
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end    
    if issparse(sce.X)
        X = single(full(sce.X));         
    else
        X = single(sce.X);
    end
    save('X.mat', '-v7.3', 'X');
    Xnorm = single(sc_norm(full(sce.X)));
    save('Xnorm.mat','-v7.3',"Xnorm","modeldir");
    g = sce.g;
    writetable(table(g),'g.csv','WriteVariableNames',false);
    sce.c_cell_id = matlab.lang.makeUniqueStrings(sce.c_cell_id);
    T = pkg.makeattributestable(sce);
    writetable(T,'c.csv');
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
    gui.gui_waitbar(fw, [], [], sprintf('Running %s...', 'py\_scimilarity'));
end

codefullpath = fullfile(codepth,'script.py');
pkg.i_addwd2script(codefullpath, wkdir, 'python');
cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
disp(cmdlinestr)
[status] = system(cmdlinestr, '-echo');
% [status2] = movefile('output.h5ad',fname);
if status == 0 && isvalid(fw)
    gui.gui_waitbar(fw, [], 'output.h5ad is written.');
end
if status == 0 && exist('output.h5ad', 'file')        
    cL = h5read('output.h5ad','/obs/predictions_unconstrained/categories');
    c = h5read('output.h5ad','/obs/predictions_unconstrained/codes');
    if any(c==0)
        cL = [cL; "undetermined"];
        c(c==0) = numel(cL);
    end
    assert(sce.NumCells==numel(c));
    if ~(isscalar(unique(sce.c_cell_type_tx)) && unique(sce.c_cell_type_tx)=="undetermined")
        sce.list_cell_attributes = [sce.list_cell_attributes, {'old_cell_type', sce.c_cell_type_tx}];
    end
    sce.c_cell_type_tx = cL(c);
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);

end


function selectedDir = selectFolder()
    % selectFolder - Prompts the user to select a folder and returns the folder path
    %
    % Output:
    %   selectedDir - The full path of the selected folder as a string. 
    %                 If the user cancels the selection, it returns an empty string.

    % Prompt title for folder selection
    promptTitle = 'Select a folder that contains the model';

    % Open a folder selection dialog box
    selectedDir = uigetdir(pwd, promptTitle);

    % Check if the user canceled the selection
    if selectedDir == 0
        fprintf('Folder selection canceled.\n');
        selectedDir = '';
    else
        fprintf('Selected folder: %s\n', selectedDir);
    end
end
