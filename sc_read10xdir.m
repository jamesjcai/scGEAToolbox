function [X,genelist,celllist,ftdone]=sc_read10xdir(selpath,coln)
%Read 10x folder
%Read files from a 10x Genomics cellranger output folder
%[X,genelist,celllist]=sc_read10xdir(selpath,coln);
%[X,genelist,celllist]=sc_read10xdir(pwd(),2);

if nargin<2, coln=2; end
if nargin<1, selpath = uigetdir; end
fprintf('Processing %s...\n',selpath);
mmfname=fullfile(selpath,'matrix.mtx');
zmmfname=fullfile(selpath,'matrix.mtx.gz');
if ~exist(mmfname,'file')
    if ~exist(zmmfname,'file')
        error('[sc_read10xdir] No matrix.mtx file.');
    else
        [~,nametxt]=fileparts(zmmfname);
        fprintf('Unzipping %s.gz...\n',nametxt);
        gunzip(zmmfname);
    end
end

ftdone=false;
ftfname=fullfile(selpath,'features.tsv');
zftfname=fullfile(selpath,'features.tsv.gz');
if ~exist(ftfname,'file')
    if ~exist(zftfname,'file')
        % error('No features.tsv file.');
    else
        [~,nametxt]=fileparts(zftfname);
        fprintf('Unzipping %s.gz...\n',nametxt);
        gunzip(zftfname);
        ftdone=true;
    end
else
    ftdone=true;
end

if ~ftdone

ftfname=fullfile(selpath,'genes.tsv');
zftfname=fullfile(selpath,'genes.tsv.gz');
if ~exist(ftfname,'file')
    if ~exist(zftfname,'file')
        % error('No features.tsv file.');
    else
        [~,nametxt]=fileparts(zftfname);
        fprintf('Unzipping %s.gz...\n',nametxt);
        gunzip(zftfname);
        ftdone=true;
    end
else
    ftdone=true;
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
    celllist=[];
end
fprintf('done.\n');
if exist(zmmfname,'file') && exist(mmfname,'file')
    delete(mmfname);
end
end
