function [T] = py_scTenifoldCko_path(sce_ori, celltype1, celltype2, targetg, ...
                                targetpathid, wkdir, ...
                                isdebug, prepare_input_only, parentfig)

T = [];
if nargin < 9, parentfig = []; end
if nargin < 8, prepare_input_only = false; end
if nargin < 7, isdebug = true; end
if nargin < 6, wkdir = []; end
if nargin < 5 || isempty(targetpathid), targetpathid = [1 3]; end
if nargin < 4 || isempty(targetg), targetg = sce_ori.g(1); end
if nargin < 3, error('Usage: [T] = py_scTenifoldCko_path(sce, celltype1, celltype2, ["ligandgene", "receptorgene"])'); end
twosided = true;
[~, targetgid]=ismember(targetg,sce_ori.g);
assert(numel(targetgid)==2)

sce = copy(sce_ori);
% -----------------
oldpth = pwd();
cleanupCwd = onCleanup(@() cd(oldpth));
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, '..', 'external', 'py_scTenifoldCko');

if isempty(wkdir) || ~isfolder(wkdir)
    wkdir = pkg.i_tempdirfile();
    cd(wkdir);
else
    disp('Using working directory provided.');
    cd(wkdir);
end

if ~prepare_input_only
 x = pyenv;
end


tmpfilelist = {'X1.mat', 'X2.mat', 'g1.txt', 'c1.txt', 'g2.txt', 'c2.txt', 'output.txt', ...
        '1/gene_name_Source.tsv', '1/gene_name_Target.tsv', ...
        '2/gene_name_Source.tsv', '2/gene_name_Target.tsv', ...
        '1/pcnet_Source.mat', '1/pcnet_Target.mat', ...
        '2/pcnet_Source.mat', '2/pcnet_Target.mat'};

pkg.i_deletefiles(tmpfilelist);   % always clear stale files, so a failed
% run cannot leave a previous run's output to be picked up as this one's

in_prepareX12intact(sce);

% One progress bar for the whole run, on the app window. It used to be
% three unparented bars in a row (one of them just a PAUSE(3)), and a
% failure left the last one open under the error dialog: the cleanup
% closes it on every exit, and errors now go to the caller, which reports
% them once.
fw = gui.myWaitbar(parentfig, [], [], 'Step 1 of 2: Building networks...');
closeFw = onCleanup(@() gui.myWaitbar(parentfig, fw, true));
in_prepareA12intact(sce);

codefullpath = fullfile(codepth,'script_path.py');
pkg.i_addwd2script(codefullpath, wkdir, 'python');

if ~prepare_input_only
    gui.myWaitbar(parentfig, fw, false, [], 'Step 2 of 2: Running scTenifoldCko.py...');
    cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
    disp(cmdlinestr)
    % https://www.mathworks.com/matlabcentral/answers/334076-why-does-externally-called-exe-using-the-system-command-freeze-on-the-third-call
    [status] = system(cmdlinestr, '-echo');
end
gui.myWaitbar(parentfig, fw);

if ~prepare_input_only

    if status == 0 && exist('output1.txt', 'file')
        T = readtable('output1.txt');
        if exist('output2.txt', 'file')
            T2 = readtable('output2.txt');
            T = {T, T2};
        end
    else
        if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
        error('scTenifoldCko runtime error.');
    end
    end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end


% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------

function in_prepareX12intact(sce)
        for id = 1:2
            if ~exist(sprintf('%d', id), 'dir')
                mkdir(sprintf('%d', id));
            end
            idx = sce.c_cell_type_tx == celltype1 | sce.c_cell_type_tx == celltype2;
            sce = sce.selectcells(idx); % OK
            sce.c_batch_id = sce.c_cell_type_tx;
            sce.c_batch_id(sce.c_cell_type_tx == celltype1) = "Source";
            sce.c_batch_id(sce.c_cell_type_tx == celltype2) = "Target";
            if issparse(sce.X)
                X = single(full(sce.X));
            else
                X = single(sce.X);
            end
            save(sprintf('X%d.mat', id), '-v7.3', 'X','targetgid');
            writematrix(sce.g, sprintf('g%d.txt', id));
            writematrix(sce.c_batch_id, sprintf('c%d.txt', id));
            fprintf('Input X%d g%d c%d written.\n', id, id, id);
            t = table(sce.g, sce.g, 'VariableNames', {' ', 'gene_name'});
            writetable(t, sprintf('%d/gene_name_Source.tsv', id), ...
                'filetype', 'text', 'Delimiter', '\t');
            writetable(t, sprintf('%d/gene_name_Target.tsv', id), ...
                'filetype', 'text', 'Delimiter', '\t');
            disp('Input gene_names written.');
        end
    end

function in_prepareA12intact(sce)
        % Log-normalised input, as in py_scTenifoldXct; Python uses these
        % networks as-is (rebuild_GRN=False).
        disp('Building A1 network...')
        X1 = sce.X(:, sce.c_cell_type_tx == celltype1);
        A = ten.i_pcnet(ten.i_lognorm(X1), 3, 0.75, false, false, symmetrize=false);
        disp('A1 network built.')
        save(sprintf('%d/pcnet_Source.mat', 1), 'A', '-v7.3');
        save(sprintf('%d/pcnet_Source.mat', 2), 'A', '-v7.3');
        disp('Building A2 network...');
        X2 = sce.X(:, sce.c_cell_type_tx == celltype2);
        A = ten.i_pcnet(ten.i_lognorm(X2), 3, 0.75, false, false, symmetrize=false);
        disp('A2 network built.');
        save(sprintf('%d/pcnet_Target.mat', 1), 'A', '-v7.3');
        save(sprintf('%d/pcnet_Target.mat', 2), 'A', '-v7.3');
    end

end
