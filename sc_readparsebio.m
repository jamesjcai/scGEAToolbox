function [X, genelist, celllist, ftdone] = sc_readparsebio(selpath, ~)

if nargin < 2, coln = 2; end
if nargin < 1, selpath = uigetdir; end
% if isempty(selpath) || selpath==0 || ~isfolder(selpath)
%     error('Need valide folder name.');
% end
fprintf('Processing %s...\n', selpath);
[out, aff] = i_guessmtxfile(selpath);
if ~isempty(out)
    disp('Found DGE.mtx');
end

if ~isempty(aff)
    mmfname = fullfile(selpath, sprintf('%sDGE.mtx', aff));
    zmmfname = fullfile(selpath, sprintf('%sDGE.mtx.gz', aff));
else
    mmfname = fullfile(selpath, out);
    zmmfname = fullfile(selpath, sprintf('%s.gz', out));
end
if ~exist(mmfname, 'file')
    if ~exist(zmmfname, 'file')
        error('No DGE.mtx file.');
    else
        [~, nametxt] = fileparts(zmmfname);
        fprintf('Unzipping %s.gz...\n', nametxt);
        gunzip(zmmfname);
    end
end

ftdone = false;
if ~isempty(aff)
    ftfname = fullfile(selpath, sprintf('%sall_genes.csv', aff));
    zftfname = fullfile(selpath, sprintf('%sall_genes.csv.gz', aff));
else
    ftfname = fullfile(selpath, 'all_genes.csv');
    zftfname = fullfile(selpath, 'all_genes.csv.gz');
end
if ~exist(ftfname, 'file')
    if ~exist(zftfname, 'file')
        % error('No features.tsv file.');
        ftdone = false;
    else
        [~, nametxt] = fileparts(zftfname);
        fprintf('Unzipping %s.gz...\n', nametxt);
        gunzip(zftfname);
        ftdone = true;
    end
else % ftfname exisiting
    ftdone = true;
end

if ~isempty(aff)
    bcfname = fullfile(selpath, sprintf('%scell_metadata.csv', aff));
    zbcfname = fullfile(selpath, sprintf('%scell_metadata.csv.gz', aff));
else
    bcfname = fullfile(selpath, 'cell_metadata.csv');
    zbcfname = fullfile(selpath, 'cell_metadata.csv.gz');
end
if ~exist(bcfname, 'file')
    if ~exist(zbcfname, 'file')
        warning('No cell_metadata.csv file.');
    else
        [~, nametxt] = fileparts(zbcfname);
        fprintf('Unzipping %s.gz...\n', nametxt);
        gunzip(zbcfname);
    end
end


if ~exist(mmfname, 'file'), error('No matrix file'); end
if ~exist(ftfname, 'file'), error('No feature file'); end

fprintf('Reading matrix file...');

% if exist(bcfname,'file')
%     [X,genelist,celllist]=sc_readmtxfile(mmfname,ftfname,bcfname,coln);
% else
%     [X,genelist]=sc_readmtxfile(mmfname,ftfname,[],coln);
%     celllist=[];
% end
[X] = sc_readmtxfile(mmfname);
if ~issparse(X)
    X = uint16(X);
end
X = X.';
fprintf('done.\n');

T1 = readtable(ftfname, 'ReadVariableNames', true, ...
    'filetype', 'text', 'Delimiter', {'\t', ',', ' ', ';', '|'}, ...
    'VariableNamingRule', 'modify');
genelist = string(T1.gene_name);

T2 = readtable(bcfname, 'ReadVariableNames', true, ...
    'filetype', 'text', 'Delimiter', {'\t', ',', ' ', ';', '|'}, ...
    'VariableNamingRule', 'modify');
celllist = string(T2.bc_wells);
ftdone = true;

if exist(zmmfname, 'file') && exist(mmfname, 'file')
    delete(mmfname);
end
if exist(zftfname, 'file') && exist(ftfname, 'file')
    delete(ftfname);
end
if exist(zbcfname, 'file') && exist(bcfname, 'file')
    delete(bcfname);
end
end

function [out, aff] = i_guessmtxfile(selpath)
out = [];
aff = [];
a = dir(selpath);
for k = 1:length(a)
    if contains(a(k).name, 'DGE.mtx')
        out = a(k).name;
        aff = extractBefore(out, 'DGE.mtx');
        continue;
    end
end
if isempty(out)
    for k = 1:length(a)
        if contains(a(k).name, 'count_matrix.mtx')
            out = 'count_matrix.mtx';
            continue;
        end
    end
end
end