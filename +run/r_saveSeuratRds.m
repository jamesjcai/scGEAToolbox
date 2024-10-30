function [status] = r_saveSeuratRds(sce, filename, wkdir)

if nargin < 3, wkdir = tempdir; end
[status] = 0;
isdebug = false;
if nargin < 2, error('run.saveSeuratRds(sce,filename)'); end
oldpth = pwd();
[isok, msg, codepath] = commoncheck_R('R_SeuratSaveRds');
if ~isok, error(msg);
    return;
end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'input.h5', 'output.Rds'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

%if ~strcmp(unique(sce.c_cell_type_tx), "undetermined")
    pkg.e_writeh5(full(sce.X), sce.g, 'input.h5', sce.c_cell_type_tx);
%else
%    pkg.e_writeh5(full(sce.X), sce.g, 'input.h5');
%end

%sc_writefile('input.txt',sce.X,sce.g);
%    if isdebug, return; end
Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end
codefullpath = fullfile(codepath,'script.R');
pkg.RunRcode(codefullpath, Rpath);

[status] = copyfile('output.Rds', filename, 'f');
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end






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
