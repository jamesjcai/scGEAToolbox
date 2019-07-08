function [X,genelist]=sc_readfile(filename,varargin)
if nargin<1
[filename, pathname] = uigetfile( ...
       {'*.csv;*.tsv;*.tab;*.txt', 'Expression Matrix Files (*.csv, *.tsv, *.tab, *.txt)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a Exprssion Matrix file');
	if isequal(filename,0), X=[]; genelist=[]; return; end
	filename=[pathname,filename];
end

if exist(filename,'file') ~= 2
    error(message('FileNotFound'));        
end

p = inputParser;
defaultType = 'tsv';
validTypes = {'tsv','mtx'};
checkType = @(x) any(validatestring(x,validTypes));

% addRequired(p,'X',@isnumeric);
addOptional(p,'type',defaultType,checkType)
parse(p,varargin{:})

switch p.Results.type
    case 'tsv'
        [X,genelist]=sc_readtsvfile(filename);
    case 'mtx'
        [X,genelist]=sc_readmtxfile(filename);
end
end