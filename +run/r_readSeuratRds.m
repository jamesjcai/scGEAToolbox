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
    X = h5read('output.h5', '/X');
elseif exist('X.csv', 'file')
    X = readmatrix('X.csv');
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
    if sce.NumCells == length(string(t.x))
        sce.c_cell_type_tx = string(t.x);
    end
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end