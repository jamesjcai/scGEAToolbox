function [X, genelist, celllist, ftdone] = sc_read10xdir(selpath, coln)
%Read 10x folder
%Read files from a 10x Genomics cellranger output folder
%[X,genelist,celllist]=sc_read10xdir(selpath,coln);
%[X,genelist,celllist]=sc_read10xdir(pwd(),2);


if nargin < 2, coln = 2; end
if nargin < 1, selpath = uigetdir; end
% if isempty(selpath) || selpath==0 || ~isfolder(selpath)
%     error('Need valide folder name.');
% end
fprintf('Processing %s...\n', selpath);
[out, aff] = i_guessmtxfile(selpath);
if isempty(out)
    selpath = fullfile(selpath, 'filtered_feature_bc_matrix');
    [out, aff] = i_guessmtxfile(selpath);
    if ~isempty(out)
        disp('Found folder ''filtered_feature_bc_matrix''');
    end
end

if ~isempty(aff)
    mmfname = fullfile(selpath, sprintf('%smatrix.mtx', aff));
    zmmfname = fullfile(selpath, sprintf('%smatrix.mtx.gz', aff));
else
    mmfname = fullfile(selpath, 'matrix.mtx');
    zmmfname = fullfile(selpath, 'matrix.mtx.gz');
end
if ~exist(mmfname, 'file')
    if ~exist(zmmfname, 'file')
        error('[sc_read10xdir] No matrix.mtx file.');
    else
        [~, nametxt] = fileparts(zmmfname);
        fprintf('Unzipping %s.gz...\n', nametxt);
        gunzip(zmmfname);
    end
end

ftdone = false;
if ~isempty(aff)
    ftfname = fullfile(selpath, sprintf('%sfeatures.tsv', aff));
    zftfname = fullfile(selpath, sprintf('%sfeatures.tsv.gz', aff));
else
    ftfname = fullfile(selpath, 'features.tsv');
    zftfname = fullfile(selpath, 'features.tsv.gz');
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

if ~ftdone
    if ~isempty(aff)
        ftfname = fullfile(selpath, sprintf('%sgenes.tsv', aff));
        zftfname = fullfile(selpath, sprintf('%sgenes.tsv.gz', aff));
    else
        ftfname = fullfile(selpath, 'genes.tsv');
        zftfname = fullfile(selpath, 'genes.tsv.gz');
    end
    if ~exist(ftfname, 'file')
        if ~exist(zftfname, 'file')
            error('No genes/features.tsv file.');
        else
            [~, nametxt] = fileparts(zftfname);
            fprintf('Unzipping %s.gz...\n', nametxt);
            gunzip(zftfname);
            ftdone = true;
        end
    else
        ftdone = true;
    end
end

if ~isempty(aff)
    bcfname = fullfile(selpath, sprintf('%sbarcodes.tsv', aff));
    zbcfname = fullfile(selpath, sprintf('%sbarcodes.tsv.gz', aff));
else
    bcfname = fullfile(selpath, 'barcodes.tsv');
    zbcfname = fullfile(selpath, 'barcodes.tsv.gz');
end
if ~exist(bcfname, 'file')
    if ~exist(zbcfname, 'file')
        warning('No barcodes.tsv file.');
    else
        [~, nametxt] = fileparts(zbcfname);
        fprintf('Unzipping %s.gz...\n', nametxt);
        gunzip(zbcfname);
    end
end


if ~exist(mmfname, 'file'), error('No matrix file'); end
if ~exist(ftfname, 'file'), error('No feature file'); end

fprintf('Reading matrix file...');
if exist(bcfname, 'file')
    [X, genelist, celllist] = sc_readmtxfile(mmfname, ftfname, bcfname, coln);
else
    [X, genelist] = sc_readmtxfile(mmfname, ftfname, [], coln);
    celllist = [];
end
fprintf('done.\n');
%{
if exist(zmmfname, 'file') && exist(mmfname, 'file')
    delete(mmfname);
end
if exist(zftfname, 'file') && exist(ftfname, 'file')
    delete(ftfname);
end
if exist(zbcfname, 'file') && exist(bcfname, 'file')
    delete(bcfname);
end
%}
end

function [out, aff] = i_guessmtxfile(selpath)
out = [];
aff = [];
a = dir(selpath);
for k = 1:length(a)
    if contains(a(k).name, 'matrix.mtx')
        out = a(k).name;
        aff = extractBefore(out, 'matrix.mtx');
        % continue;
        break;
    end
end
end