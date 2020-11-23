function [X,genelist,celllist]=sc_read10xdir(selpath,coln)
%Read files from a 10x Genomics cellranger output folder

if nargin<2, coln=2; end
if nargin<1, selpath = uigetdir; end
fprintf('Processing %s...\n',selpath);
mmfname=fullfile(selpath,'matrix.mtx');
zmmfname=fullfile(selpath,'matrix.mtx.gz');
if ~exist(mmfname,'file')
    if ~exist(zmmfname,'file')
        error('No matrix file.');
    else
        [~,nametxt]=fileparts(zmmfname);
        fprintf('Unzipping %s.gz...\n',nametxt);
        gunzip(zmmfname);
    end
end
ftfname=fullfile(selpath,'features.tsv');
zftfname=fullfile(selpath,'features.tsv.gz');
if ~exist(ftfname,'file')
    if ~exist(zftfname,'file')
        error('No features.tsv file.');
    else
        [~,nametxt]=fileparts(zftfname);
        fprintf('Unzipping %s.gz...\n',nametxt);
        gunzip(zftfname);
    end
end

bcfname=fullfile(selpath,'barcodes.tsv');
zbcfname=fullfile(selpath,'barcodes.tsv.gz');
if ~exist(bcfname,'file')
    if ~exist(zbcfname,'file')
        warning('No barcodes.tsv file.');
    else
        [~,nametxt]=fileparts(zbcfname);
        fprintf('Unzipping %s.gz...\n',nametxt);        
        gunzip(zbcfname);
    end
end




if ~exist(mmfname,'file')
    error('No matrix file');
end
if ~exist(ftfname,'file')
    error('No feature file');
end

fprintf('Reading matrix file...');
if exist(bcfname,'file')    
    [X,genelist,celllist]=sc_readmtxfile(mmfname,ftfname,bcfname,coln);
else
    [X,genelist]=sc_readmtxfile(mmfname,ftfname,[],coln);
end
fprintf('done.\n');
if exist(zmmfname,'file') && exist(mmfname,'file')
    delete(mmfname);
end
end
