function [status] = r_saveSeuratRds(sce, filename, wkdir)

if nargin < 3, wkdir = pkg.i_tempdirfile(); end
[status] = 0;
isdebug = false;

if nargin < 2, error('run.r_saveSeuratRds(sce,filename)'); end
oldpth = pwd();
cleanupCwd = onCleanup(@() cd(oldpth));
[isok, msg, codepath] = commoncheck_R('R_SeuratSaveRds');
if ~isok, error('%s', msg); end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'input.h5', 'output.Rds'};
pkg.i_deletefiles(tmpfilelist);   % always clear stale files, so a failed
% run cannot leave a previous run's output to be picked up as this one's

% if ~strcmp(unique(sce.c_cell_type_tx), "undetermined")
% Cell IDs become the Seurat colnames; script.R falls back to C1..Cn without them
cellid = string(sce.c_cell_id(:));
if numel(cellid) ~= sce.NumCells || any(ismissing(cellid) | cellid == "")
    cellid = [];
end
pkg.e_writeh5(full(sce.X), sce.g, 'input.h5', sce.c_cell_type_tx, sce.c_batch_id, cellid);

Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end
codefullpath = fullfile(codepath,'script.R');
pkg.i_runrcode(codefullpath, Rpath);

% output.Rds was deleted above, so its absence means R failed; COPYFILE
% would only return 0, which the caller used to ignore.
if ~isfile('output.Rds')
    error('run:r_saveSeuratRds:noOutput', ...
        'R did not write the Seurat file. The R output in the Command Window should say why.');
end
[status] = copyfile('output.Rds', filename, 'f');
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
end


%{
    function [indptr, indices, data] = convert_sparse_to_indptr(X)

        % Check if X is sparse
        if ~issparse(X)
            error('Input matrix X must be a sparse matrix');
        end

        % Get matrix dimensions
        [~, n] = size(X);

        % Initialize indptr with 1 and n+1
        indptr = [1, n+1];

        % Find non-zero elements and their indices
        [row, col] = find(X);

        % Sort by columns for efficient construction
        [~, sort_idx] = sort(col);
        row = row(sort_idx);
        col = col(sort_idx);

        % Accumulate column counts for indptr
        for i = 1:n
            indptr(i+1) = indptr(i) + sum(col == i);
        end

        % Assign indices and data
        indices = row;
        data = full(X(row, col));  % Extract non-zero values

    end
%}
