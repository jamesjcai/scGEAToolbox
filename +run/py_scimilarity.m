function [celltypes, stats] = py_scimilarity(sce, modeldir, wkdir, ...
    target_celltypes, isdebug, prepare_input_only, cellidx)
%PY_SCIMILARITY Annotate cells against a SCimilarity reference model.
%
%   celltypes = run.py_scimilarity(sce, modeldir, wkdir)
%   celltypes = run.py_scimilarity(sce, modeldir, wkdir, target_celltypes, ...
%       isdebug, prepare_input_only, cellidx)
%
%   CELLIDX restricts the run to those columns of SCE.X, and CELLTYPES
%   then carries one label per index, in the order given. Leave it empty
%   for every cell.
%
%   STATS is one row per returned label, with the confidence columns
%   get_predictions_knn produces and this function used to throw away:
%   vsAll is the winning label's share of the 50 reference neighbours,
%   vs2nd its share against the runner-up alone, and the _weighted
%   variants are those counts under the 1/distance weighting the
%   prediction itself uses. min_dist and max_dist bound the distance to
%   those neighbours. STATS is empty when the working folder holds an
%   older script_mat.py that does not write it.
%
%   Subsetting here rather than in the caller is deliberate. SCimilarity
%   embeds and queries each cell independently, and normalisation is per
%   cell too, so a subset run returns exactly the labels a full run would
%   have given those cells. It also holds the FULL() below to the cells
%   actually being sent, which on a large dataset is the difference
%   between a dense copy that fits in memory and one that does not.
%
%   Example:
%     cell_type = run.py_scimilarity(sce, 'Y:\jcai\models\model_v1.1', ...
%         'C:\Users\jcai\Downloads');

celltypes = [];
stats = table();
if nargin < 7, cellidx = []; end
if nargin < 6, prepare_input_only = false; end
if nargin < 5, isdebug = true; end
if nargin < 4, target_celltypes = ''; end
if nargin < 3, wkdir = pkg.i_tempdirfile(); end
if nargin < 2, modeldir = selectFolder; end
if isempty(modeldir) || ~exist(modeldir, 'dir')
    error('Model folder does not exist or is invalid.');
end

oldpth = pwd();
cleanupCwd = onCleanup(@() cd(oldpth));
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, '..', 'external', 'py_scimilarity');

if isempty(wkdir) || ~isfolder(wkdir)
    cd(codepth);
else
    disp('Using working directory provided.');
    cd(wkdir);
end
% winopen(wkdir);

% fw = gui.gui_waitbar([], [], 'Checking Python environment...');

x = pyenv;
try
    pkg.i_add_conda_python_path;
catch
    % best-effort: fall back to default pyenv if conda path not found
end
codepth = pkg.i_normalizepath(codepth);


if ~prepare_input_only
        codefullpath = fullfile(codepth,'require.py');
        cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);

        disp(cmdlinestr)
        [status, cmdout] = system(cmdlinestr, '-echo');
        if status ~= 0
            error('%s', cmdout);
        else
            disp('Code requirement check is done.')
        end
    end

% try
pkg.i_deletefiles({'input.h5ad', 'output.h5ad','tg.csv'});
tmpfilelist = {'Xnorm.mat', 'X.mat', 'g.csv', 'c.csv', 'tg.csv', ...
        'input.h5ad', 'output.h5ad', 'output.csv', 'output_stats.csv'};
pkg.i_deletefiles(tmpfilelist);   % always clear stale files, so a failed
% run cannot leave a previous run's output to be picked up as this one's

% Take the requested cells first. The two steps below are per gene and per
% cell respectively, so they commute with a column subset, and doing it
% here keeps the gene collapse and the FULL() off cells nobody asked for.
X = sce.X;
if ~isempty(cellidx)
    X = X(:, cellidx);
end

% align_dataset builds the aligned matrix with anndata.concat, which refuses
% a duplicated var index. The caller upper-cases the gene list to match
% SCimilarity's upper-case HGNC gene_order, and that alone can put two genes
% on one symbol (Pisd and PISD). Collapse repeats by summing their raw
% counts, which is the usual merge and leaves every cell's library size
% unchanged, so the normalisation in script_mat.py is unaffected.
[genelist, X] = i_collapseduplicategenes(sce.g, X);

% RAW COUNTS GO TO PYTHON, not normalised values. This function used to
% library-size normalise here, which put the normalisation BEFORE
% align_dataset drops every gene outside the model's 28k gene space - so the
% library size counted genes the model never sees, and every value the
% encoder received was off by a data-dependent factor. Both SCimilarity
% tutorials align first and call lognorm_counts second, and script_mat.py now
% does the same. On the repo's example data the genes outside the model space
% carry 1.8-2.0% of counts, but new_example_sce.mat has one cell at 83%.
Xcounts = single(full(X));

    % if ~isempty(target_celltypes)
    %    writetable(table(target_celltypes),'tg.csv','WriteVariableNames',false);
    % end
    % g = sce.g;
    % writetable(table(g),'g.csv','WriteVariableNames',false);
    % sce.c_cell_id = matlab.lang.makeUniqueStrings(sce.c_cell_id);
    % T = pkg.i_makeattributestable(sce);
    % writetable(T,'c.csv');


g = cellstr(genelist);
if ~isempty(target_celltypes)
        tg = cellstr(target_celltypes);
        save('X.mat','-v7.3',"Xcounts","modeldir","tg","g");
    else
        save('X.mat','-v7.3',"Xcounts","modeldir","g");
    end
% catch ME
%     if pkg.i_isvalid(fw)
%          gui.gui_waitbar(fw, true);
%     end
%     errordlg(ME.message,'');
%     return;
% end
% if pkg.i_isvalid(fw)
%     gui.gui_waitbar(fw, [], [], 'Checking Python environment is complete');
%     pause(0.5);
%     gui.gui_waitbar(fw, [], [], sprintf('Running %s...', 'py\_scimilarity'));
% end
codefullpath = fullfile(codepth,'script_mat.py');
pkg.i_addwd2script(codefullpath, wkdir, 'python');
% The copy i_addwd2script leaves in wkdir is what the user runs by hand after
% a prepare-input-only run, and it imports _scimilarity_env by name. Put that
% next to it so the prepared folder stands on its own.
copyfile(fullfile(codepth, '_scimilarity_env.py'), wkdir);
cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
disp(cmdlinestr)

if ~prepare_input_only
    [status] = system(cmdlinestr, '-echo');
    % [status2] = movefile('output.h5ad',fname);
    % if status == 0 && pkg.i_isvalid(fw)
    %     gui.gui_waitbar(fw, [], 'output.csv is written.');
    % end
    if status == 0 && exist('output.csv', 'file')
        t = readtable('output.csv','ReadVariableNames', true, ...
            'VariableNamingRule', 'modify');
        celltypes = string(t.x0);
        % Confidence per label, which get_predictions_knn computes anyway
        % and this function used to discard. Guarded rather than assumed:
        % a working folder prepared by an older release still has a
        % script_mat.py that writes no stats file.
        if exist('output_stats.csv', 'file')
            stats = readtable('output_stats.csv', 'ReadVariableNames', true, ...
                'VariableNamingRule', 'modify');
            stats = removevars(stats, 1);   % the pandas index column
        end
        % cL = h5read('output.h5ad','/obs/predictions_unconstrained/categories');
        % c = h5read('output.h5ad','/obs/predictions_unconstrained/codes');
        % if any(c==0)
        %     cL = [cL; "undetermined"];
        %     c(c==0) = numel(cL);
        % end
    end
else
    % script_mat.py, not script.py: I_ADDWD2SCRIPT copies script_mat.py into
    % the working folder, and it is the one that reads the X.mat written
    % above. script.py is a stale variant expecting c.csv and g.csv, which
    % nothing writes any more.
    disp('Input files are prepared. To do the analysis, run script_mat.py in the working folder.')
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

end


function [genelist, X] = i_collapseduplicategenes(genelist, X)
% Sum the rows of X whose gene symbols repeat, keeping first-seen order.

genelist = string(genelist);
[genelist, ~, idx] = unique(genelist, 'stable');
nduplicate = numel(idx) - numel(genelist);
if nduplicate == 0, return; end

fprintf('Collapsing %d duplicate gene symbol(s) by summing counts.\n', ...
    nduplicate);
selector = sparse(idx, (1:numel(idx))', 1, numel(genelist), numel(idx));
X = selector*X;
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
