function [sce, metadata] = r_readSeuratRds(filename, wkdir)

if nargin < 2, wkdir = pkg.i_tempdirfile(); end
sce = [];
metadata = [];

if nargin < 1, error('run.r_readSeuratRds(filename)'); end
oldpth = pwd();
cleanupCwd = onCleanup(@() cd(oldpth));
[isok, msg, codepth] = commoncheck_R('R_SeuratReadRds');
if ~isok, error('%s', msg); end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

isdebug = false;
tmpfilelist = {'inputrdsfile.txt', 'output.h5', ...
'g.csv', 'X.csv', 'umap.csv', 'barcodes.csv', ...
'annotation.csv', 'metadata.csv'};
pkg.i_deletefiles(tmpfilelist);   % always clear stale files, so a failed
% run cannot leave a previous run's output to be picked up as this one's


writematrix(filename, 'inputrdsfile.txt');
Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end
codefullpath = fullfile(codepth,'script.R');
pkg.i_runrcode(codefullpath, Rpath);

g = [];
if exist('g.csv', 'file')
    t = readtable('g.csv', 'Delimiter', ',', 'VariableNamingRule', 'modify');
    if isa(t, 'table')
        g = string(t.x);
    end
end

if exist('output.h5', 'file')
    % X = h5read('output.h5', '/X');
    filenm = 'output.h5';
    grouptag = '/';
    data = pkg.e_guessh5field(filenm, {grouptag}, {'data'}, true);
    indices = pkg.e_guessh5field(filenm, {grouptag}, {'indices'}, true);
    indptr = pkg.e_guessh5field(filenm, {grouptag}, {'indptr'}, true);
    shape = pkg.e_guessh5field(filenm, {grouptag}, {'shape'}, true);


    if ~isMATLABReleaseOlderThan('R2025a')
        X = spalloc(shape(1), shape(2), length(data), 'single');
    else
        X = spalloc(shape(1), shape(2), length(data));
    end

    % X = spalloc(shape(1), shape(2), length(data));

    for k = 1:length(indptr) - 1
        i = indptr(k) + 1:indptr(k+1);
        y = indices(i) + 1;
        X(y, k) = data(i);
    end
    if isequal(size(X), shape)
        warning('Matrix size changed.');
    end
    disp('X is read.');
else
    error('run.r_readSeuratRds:noOutput', ...
        ['R finished but did not write %s to %s. The R console ', ...
         'output above should say why.'], 'output.h5', pwd);
end

X = pkg.e_uint2sparse(X);
sce = SingleCellExperiment(X, g);

barcodesfile = 'barcodes.csv';
if exist(barcodesfile, 'file')
    disp('Reading c_cell_id from barcodes.csv');
    t = readtable('barcodes.csv', 'Delimiter', ',', 'VariableNamingRule', 'modify');
    if ~isempty(t) && ismember('x', t.Properties.VariableNames)
        id = string(t.x);
        if length(id) == sce.NumCells, sce.c_cell_id = id; end
    end
end
umapfile = 'umap.csv';
if exist(umapfile, 'file')
    disp('Reading s from umap.csv');
    t = readtable('umap.csv', 'Delimiter', ',', 'VariableNamingRule', 'modify');
    [y, idx] = ismember(string(t.Var1), sce.c_cell_id);
    if all(y)
        s = table2array(t(:, 2:end));
        sce.s = s(idx, :);
        sce.struct_cell_embeddings.umap2d = s;
    end
end
batchfile = 'batch.csv';
if exist(batchfile,'file')
    disp('Reading batchid from batch.csv');
    t=readtable('batch.csv','VariableNamingRule', 'modify');
    id=string(t.x);
    if length(id) == sce.NumCells, sce.c_batch_id = id; end
end


% if exist('annotation.csv', 'file')
%     disp('Reading celltype from annotation.csv');
%     t = readtable('annotation.csv', 'Delimiter', ',');
%     if ~isempty(t) && ismember('x', text.Properties.VariableNames)
%         if sce.NumCells == length(string(t.x))
%             sce.c_cell_type_tx = string(t.x);
%         end
%     end
% else

% if isdeployed || ~isempty(which('celltype.csv'))
celltypefile = 'celltype.csv';
if exist(celltypefile, 'file')
    disp('Reading celltype from celltype.csv');
    t = readtable('celltype.csv', 'Delimiter', ',','VariableNamingRule', 'modify');
    if ~isempty(t) && ismember('x', t.Properties.VariableNames)
        if sce.NumCells == length(string(t.x))
            sce.c_cell_type_tx = string(t.x);
        end
    end
end
metadatafile = 'metadata.csv';
if exist(metadatafile, 'file')
    try
    disp('Reading metadata from metadata.csv');
        t = readtable('metadata.csv', 'Delimiter', ',', ...
            'VariableNamingRule', 'modify');
        if sce.NumCells == height(t)
            metadata = in_alignbybarcode(t, sce.c_cell_id);
        end
    catch
        % metadata is optional; sce is still valid without it
    end
end

% Keep the per-cell metadata with the SCE, not just in the caller's workspace,
% so it stays viewable and usable as grouping variables after import.
if ~isempty(metadata)
    sce = in_addmetadataattributes(sce, metadata);
end

% The R side only writes celltype.csv for a column literally named CellType or
% celltype. When that misses, infer the annotation column from the metadata
% table and take it only when the evidence is strong; callers with a UI can
% offer gui.i_pickcelltypecolumn for the ambiguous cases.
if isempty(sce.c_cell_type_tx) && ~isempty(metadata)
    [colname, colscore] = pkg.i_guesscelltypecol(metadata, sce.NumCells);
    if colscore >= 70
        ctype = strtrim(string(metadata.(colname)));
        ctype(ismissing(ctype) | strlength(ctype) == 0) = "undetermined";
        sce.c_cell_type_tx = ctype;
        fprintf('Cell type read from metadata column ''%s'' (score %.0f).\n', ...
            colname, colscore);
    end
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
end

function sce = in_addmetadataattributes(sce, t)
% Store each per-cell metadata column in list_cell_attributes, which makes them
% appear in the cell attribute table and in the grouping-variable selectors.
% The row name column is skipped because it is already c_cell_id, and constant
% columns are skipped because they carry no information.

varnames = string(t.Properties.VariableNames);
existing = string(sce.list_cell_attributes(1:2:end));
keep = cell(1, 2*numel(varnames));
nkeep = 0;

for k = 1:numel(varnames)
    v = t.(varnames(k));
    if size(v, 2) ~= 1, continue; end
    if any(existing == varnames(k)), continue; end
    if numel(unique(string(v))) < 2, continue; end
    if k == 1 && numel(sce.c_cell_id) == height(t) && ...
            isequal(string(v), string(sce.c_cell_id))
        continue;
    end
    keep{2*nkeep+1} = char(varnames(k));
    keep{2*nkeep+2} = v;
    nkeep = nkeep + 1;
end

if nkeep > 0
    sce.list_cell_attributes = [sce.list_cell_attributes, keep(1:2*nkeep)];
end
end

function t = in_alignbybarcode(t, cellid)
% Reorder metadata rows to match sce.c_cell_id. write.csv puts the meta.data row
% names in the first column; when those are the barcodes the join is exact,
% otherwise the original row order is kept.

if isempty(cellid) || width(t) < 1, return; end
rowid = string(t.(t.Properties.VariableNames{1}));
[isfound, idx] = ismember(string(cellid), rowid);
if all(isfound) && numel(unique(idx)) == numel(idx)
    t = t(idx, :);
end
end
