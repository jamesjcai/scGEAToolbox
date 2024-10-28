function [status] = r_copykat(sce, wkdir)

if nargin < 2
    wkdir = tempdir; 
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
    if isempty(wkdir), return; end
end

extprogname = 'R_copykat';

[status] = 0;
isdebug = true;

oldpth = pwd();
[isok, msg, codepath] = commoncheck_R('R_copykat');
if ~isok, error(msg);
    return;
end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'input.h5'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
pkg.e_writeh5(full(sce.X), sce.g, 'input.h5');
Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end
codefullpath = fullfile(codepath,'script.R');
pkg.RunRcode(codefullpath, Rpath);

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
