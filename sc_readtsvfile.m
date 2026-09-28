function [X, genelist, celllist] = sc_readtsvfile(filename, genecolnum)
% Read TSV/TXT file
if nargin < 2
    genecolnum = 1;
end
if nargin < 1
    [filename, pathname] = uigetfile( ...
        {'*.csv;*.tsv;*.tab;*.txt', ...
        'Expression Matrix Files (*.csv, *.tsv, *.tab, *.txt)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Pick a Exprssion Matrix file');
    if isequal(filename, 0), X = [];
        genelist = [];
        return;
    end
    filename = fullfile(pathname, filename);
end

if exist(filename, 'file') ~= 2
    error('FileNotFound');
end
if nargout > 2
    [X, genelist, celllist] = i_read_exprmat(filename, genecolnum);
    if size(X, 2) ~= length(celllist)
        celllist = strrep(celllist, '"', '');
        celllist(strlength(celllist) == 0) = [];
    end
    if size(X, 2) ~= length(celllist)
        warning('LENGTH(BARCODELIST) is not equal to SIZE(X,2)')
        if size(X, 2) - length(celllist) == -1
            celllist = celllist(2:end);
        end
    end
else
    [X, genelist] = i_read_exprmat(filename, genecolnum);
end

end

function [X, genelist, sampleid] = i_read_exprmat(filename, genecolnum, verbose)
% Validate input args
narginchk(1, Inf);

% Get Filename
if ~ischar(filename) && ~(isstring(filename) && isscalar(filename))
    error('FileNameMustBeString');
end
filename = char(filename);

% Make sure file exists
if exist(filename, 'file') ~= 2
    error('FileNotFound');
end


if nargin < 2, genecolnum = 1; end
if nargin < 3, verbose = true; end
if verbose, fprintf('Reading %s ...... ', filename); end
T = readtable(filename, 'filetype', 'text', 'HeaderLines', 1, ...
    'VariableNamingRule', 'modify');
X = table2array(T(:, (1 + genecolnum):end));
X = sparse(X);
if ~isMATLABReleaseOlderThan('R2025a')
    X = single(X);
end
genelist = string(table2array(T(:, 1:genecolnum)));
if nargout > 2
    fid = fopen(filename);
    sampleid = string(strsplit(fgetl(fid), {'\t', ','})');
    fclose(fid);
end
if verbose, fprintf('done.\n'); end
end