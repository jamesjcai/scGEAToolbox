function [sce] = sc_readgeoaccess(acc)

% if length(strsplit(acc,{',',';',' '}))>1
% end

url = sprintf('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s', acc);
a = webread(url);
b = strsplit(a, '\n');
c = string(b(contains(b, acc)))';
c = c(contains(c, 'ftp'));

if ~(isscalar(c) || length(c) >= 3)
    disp(url)
    % warning('Check GEO supplementary file list.');
    c = c(1);
end

barcodes = [];

if length(c) >= 3
    %switch length(c)
    %    case 3
    c1 = c(contains(c, 'mtx'));
    if isempty(c1), error('MTX file not found.'); end
    f1 = i_setupfile(c1);
    if isempty(f1), error('MTX file name not processed.'); end

    c2 = c(contains(c, 'genes'));
    if isempty(c2), c2 = c(contains(c, 'features')); end
    if isempty(c2), error('GENES/FEATURES file not found.'); end
    f2 = i_setupfile(c2);
    if isempty(f2), error('GENES/FEATURES file name not processed.'); end

    c3 = c(contains(c, 'barcodes'));
    f3 = [];
    if isempty(c3)
        warning('BARCODES file not found.');
    else
        f3 = i_setupfile(c3);
    end
    if isempty(f2)
        [X, g] = sc_readmtxfile(f1, f2);
    else
        [X, g, barcodes] = sc_readmtxfile(f1, f2, f3);
    end
elseif isscalar(c)
    txtnotfound = false;
    c1 = c(contains(c, 'txt'));
    if isempty(c1)
        c1 = c(contains(c, 'csv'));
        if isempty(c1)
            c1 = c(contains(c, 'tsv'));
            if isempty(c1)
                txtnotfound = true;
                % error('TXT/CSV/TSV file not found.');
            end
        end
    end
    if ~txtnotfound
        disp("Found TXT/CSV/TSV file.");
        f1 = i_setupfile(c1);
        if isempty(f1), error('TXT/CSV/TSV file name not processed.'); end
        [X, g] = sc_readtsvfile(f1);
    else
        c1 = c(contains(c, 'h5'));
        if isempty(c1)
            error('File not found.');
        end
        disp("Found H5 file.");
        
        f1 = i_setupfile2(c1);
        if isempty(f1), sce=[]; return; end

        if strcmpi(f1(end-2:end), '.gz')
            files=gunzip(f1,tempdir);
            [X, g, barcodes] = sc_read10xh5file(files{1});
        elseif strcmpi(f1(end-2:end), '.h5')
            [X, g, barcodes] = sc_read10xh5file(f1);
        end
    end
end


% if sum(upper(extractBefore(g,8))=='GRCH38_')>10
%     warning('Gene names contain prefix GRCH38_.');
%     try
%         g=extractAfter(g,8);
%     catch ME
%         warning('Gene names contain prefixes.');
%         rethrow(ME);
%     end
% end


if length(g) == size(X, 1)
    sce = SingleCellExperiment(X, g);
elseif length(g) == size(X, 2)
    disp('X is transposed.');
    sce = SingleCellExperiment(X.', g);
else
    warning('size(sce.X,1)~=length(sce.g)');
    sce = SingleCellExperiment(X, g);
end

metainfo = sprintf("Source: %s", acc);
sce = sce.appendmetainfo(metainfo);
fprintf(['The data was downloaded from the National Center', ...
    ' for Biotechnology Information Gene Expression Omnibus (GEO) ', ...
        'with the accession ID %s.\n'], acc);

    % function i_tryh5(c)
    %     c1=c(contains(c,'tsv'));
    % end
    if ~isempty(barcodes)
        sce.c_cell_id = barcodes;
    end
end




function f = i_setupfile(c)

% https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html#FTP
    try
        tmpd = tempdir;
        [x] = regexp(c(1), '<a href="ftp://(.*)">(ftp', 'match');
        x = string(textscan(x, '<a href="ftp://%s'));
        x = append("https://", extractBefore(x, strlength(x)-5));
        if ~(ismcc || isdeployed)
            %#exclude urldecode
            x = urldecode(x);
        else
            x = pkg.urldecoding(x);
        end
        fprintf('Downloading %s\n', x)
        files = gunzip(x, tmpd);
        f = files{1};
        if strcmpi(f(end-3:end), '.tar')
            f = untar(f, tmpd);
            if length(f) > 1
                f = [];
            end
        end
    catch ME
        disp(ME.message)
        f = [];
    end
end


function f = i_setupfile2(c)
    try
        tmpd = tempname;
        [x] = regexp(c(1), '<a href="ftp://(.*)">(ftp', 'match');
        x = string(textscan(x, '<a href="ftp://%s'));
        x = append("https://", extractBefore(x, strlength(x)-5));
        if ~(ismcc || isdeployed)
            %#exclude urldecode
            x = urldecode(x);
        else
            x = pkg.urldecoding(x);
        end
        fprintf('Downloading %s\n', x)
        f = websave(tmpd, x);
    catch
        f = [];
    end
end
