function [X, genelist] = sc_readfile(filename, varargin)
    % SC_READFILE  Load single‑cell expression data files
    %
    %   [X, genelist] = sc_readfile(filename)
    %   [X, genelist] = sc_readfile(filename, 'type', 'h5');
    %
    %   Inputs:
    %     filename — path to data file, or empty to prompt UI
    %     'type'   — 'tsv' (default), 'mtx', or 'h5'
    %
    %   Outputs:
    %     X        — expression matrix (genes × cells)
    %     genelist — cell array or string vector of gene names
    
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
    filename = [pathname, filename];
end

if exist(filename, 'file') ~= 2
    error('FileNotFound');
end

p = inputParser;
defaultType = 'tsv';
validTypes = {'tsv', 'mtx', 'h5'};
checkType = @(x) any(validatestring(x, validTypes));

% addRequired(p,'X',@isnumeric);
addOptional(p, 'type', defaultType, checkType)
parse(p, varargin{:})

switch p.Results.type
    case 'tsv'
        [X, genelist] = sc_readtsvfile(filename);
    case 'mtx'
        featurestxtfile = 'features.tsv';
        if exist(featurestxtfile, 'file')
            [X, genelist] = sc_readmtxfile(filename, featurestxtfile);
        else
            error('Need file ''features.tsv''');
        end
    case 'h5'
        [X, genelist] = sc_read10xh5file(filename);
end
end
