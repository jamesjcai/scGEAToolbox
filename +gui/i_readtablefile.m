function [t, fname] = i_readtablefile(parentfig, nrows, dlgtitle)
% I_READTABLEFILE - Pick a delimited text file and read it into a table.
%
%   [t, fname] = gui.i_readtablefile(parentfig)
%   [t, fname] = gui.i_readtablefile(parentfig, nrows, dlgtitle)
%
% Inputs:
%   parentfig - parent figure, raised again after the file dialog (default: [])
%   nrows     - required number of rows, e.g. sce.NumCells; [] accepts any
%               number of rows (default: [])
%   dlgtitle  - title of the file dialog (default: 'Pick a per-cell table file')
%
% Outputs:
%   t     - the table, [] when cancelled or when the file could not be read
%   fname - name of the file read, '' when cancelled
%
% Reading is done with import options rather than a bare readtable call for two
% reasons. readtable auto-detects only its own short list of extensions, so a
% .tsv file errors out with UnrecognizedExtension before it is ever read. And on
% a tab-delimited file readtable detects tab AND space as delimiters, which
% splits a label like "T cell" across two columns and leaves the header
% unusable; when a tab is among the delimiters found, the file is read again
% with the tab as the only one.
%
% see also: gui.i_pickworkspacetable, gui.callback_AssignCellTypeFromAttrib,
%           gui.sc_cellattribeditor

if nargin < 1, parentfig = []; end
if nargin < 2, nrows = []; end
if nargin < 3 || isempty(dlgtitle), dlgtitle = 'Pick a per-cell table file'; end

t = [];
fname = '';

[f, pathname] = uigetfile( ...
    {'*.csv;*.txt;*.tsv', 'Table Files (*.csv, *.txt, *.tsv)'; ...
    '*.*', 'All Files (*.*)'}, dlgtitle);
if pkg.i_isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure')
    figure(parentfig);
end
if isequal(f, 0), return; end

try
    t = in_readdelimited(fullfile(pathname, f));
catch ME
    gui.myErrordlg(parentfig, ME.message, ME.identifier);
    t = [];
    return;
end

if ~isempty(nrows) && height(t) ~= nrows
    gui.myErrordlg(parentfig, sprintf( ...
        'The table has %d rows; %d cells are expected.', height(t), nrows));
    t = [];
    return;
end

fname = f;
end

function t = in_readdelimited(filename)
opts = detectImportOptions(filename, 'FileType', 'text');
istab = strcmp(opts.Delimiter, '\t') | strcmp(opts.Delimiter, sprintf('\t'));
if numel(opts.Delimiter) > 1 && any(istab)
    opts = detectImportOptions(filename, 'FileType', 'text', 'Delimiter', '\t');
end
opts.VariableNamingRule = 'modify';

% Modifying headers into valid identifiers is what is being asked for here, so
% the warning readtable raises about having done it is noise. The original
% headers stay available in the table's VariableDescriptions property.
warnstate = warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
cleanupObj = onCleanup(@() warning(warnstate));

t = readtable(filename, opts);
end
