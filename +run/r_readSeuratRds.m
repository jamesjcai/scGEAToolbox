function [sce] = r_readSeuratRds(filename, wkdir)

if nargin < 2, wkdir = tempdir; end
sce = [];
if nargin < 1, error('run.r_readSeuratRds(filename)'); end
oldpth = pwd();
[isok, msg, codepth] = commoncheck_R('R_SeuratReadRds');
if ~isok, error(msg);
    sce = [];
    return;
end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

isdebug = false;
tmpfilelist = {'inputrdsfile.txt', 'output.h5', ...
    'g.csv', 'X.csv', 'umap.csv', 'barcodes.csv', 'annotation.csv'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

writematrix(filename, 'inputrdsfile.txt');
Rpath = getpref('scgeatoolbox', 'rexecutablepath');
codefullpath = fullfile(codepth,'script.R');
pkg.RunRcode(codefullpath, Rpath);

g = [];
if exist('g.csv', 'file')
    t = readtable('g.csv', 'Delimiter', ',');
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
end

X = pkg.e_uint2sparse(X);
sce = SingleCellExperiment(X, g);


if exist('barcodes.csv', 'file') && ~isempty(sce)
    t = readtable('barcodes.csv', 'Delimiter', ',');
    id = string(t.x);
    sce.c_cell_id = id;
end

if exist('umap.csv', 'file') && ~isempty(sce) && ~isempty(sce.c_cell_id)   
    t = readtable('umap.csv', 'Delimiter', ',');
    [y, idx] = ismember(string(t.Var1), sce.c_cell_id);
    if all(y)
        s = table2array(t(:, 2:end));
        sce.s = s(idx, :);
    end
end
%     if exist('batchid.csv','file') && ~isempty(sce)
%         t=readtable('batchid.csv');
%         id=string(t.x);
%         sce.c_batch_id=id;
%     end


if exist('annotation.csv', 'file') && ~isempty(sce) && ~isempty(sce.c_cell_id)
    t = readtable('annotation.csv', 'Delimiter', ',');
    if ~isempty(t) && contains(t.Properties.VariableNames,'x')
        if sce.NumCells == length(string(t.x))
            sce.c_cell_type_tx = string(t.x);
        end
    end
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end