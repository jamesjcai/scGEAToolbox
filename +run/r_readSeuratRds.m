function [sce, metadata] = r_readSeuratRds(filename, wkdir)

if nargin < 2, wkdir = tempdir; end
sce = [];
metadata = [];

if nargin < 1, error('run.r_readSeuratRds(filename)'); end
oldpth = pwd();
[isok, msg, codepth] = commoncheck_R('R_SeuratReadRds');
if ~isok, error(msg), return; end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

isdebug = false;
tmpfilelist = {'inputrdsfile.txt', 'output.h5', ...
    'g.csv', 'X.csv', 'umap.csv', 'barcodes.csv', ...
    'annotation.csv', 'metadata.csv'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end


writematrix(filename, 'inputrdsfile.txt');
Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end
codefullpath = fullfile(codepth,'script.R');
pkg.RunRcode(codefullpath, Rpath);

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
    X = spalloc(shape(1), shape(2), length(data));    
    for k = 1:length(indptr) - 1
        i = indptr(k) + 1:indptr(k+1);
        y = indices(i) + 1;
        X(y, k) = data(i);
    end
    if isequal(size(X), shape)
        warning('Matrix size changed.');
    end
    disp('X is read.');
end

X = pkg.e_uint2sparse(X);
sce = SingleCellExperiment(X, g);


if exist('barcodes.csv', 'file')
    disp('Reading c_cell_id from barcodes.csv');
    t = readtable('barcodes.csv', 'Delimiter', ',', 'VariableNamingRule', 'modify');
    if ~isempty(t) && ismember('x', t.Properties.VariableNames)
        id = string(t.x);
        if length(id) == sce.NumCells, sce.c_cell_id = id; end
    end
end

if exist('umap.csv', 'file')
    disp('Reading s from umap.csv');
    t = readtable('umap.csv', 'Delimiter', ',', 'VariableNamingRule', 'modify');
    [y, idx] = ismember(string(t.Var1), sce.c_cell_id);
    if all(y)
        s = table2array(t(:, 2:end));
        sce.s = s(idx, :);
    end
end

if exist('batch.csv','file')
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
if exist('celltype.csv', 'file')
    disp('Reading celltype from celltype.csv');
    t = readtable('celltype.csv', 'Delimiter', ',','VariableNamingRule', 'modify');
    if ~isempty(t) && ismember('x', t.Properties.VariableNames)
        if sce.NumCells == length(string(t.x))
            sce.c_cell_type_tx = string(t.x);
        end
    end    
end

if exist('metadata.csv', 'file')
    try
    disp('Reading metadata from metadata.csv');
        t = readtable('metadata.csv', 'Delimiter', ',', ...
            'VariableNamingRule', 'modify');
        if sce.NumCells == height(t)
            metadata = t;
        end
    catch
    end
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end