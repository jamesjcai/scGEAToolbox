function [X,genelist,celllist]=sc_readtsvfile(filename)
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
%     try
if nargout>2
    [X,genelist,celllist]=i_read_exprmat(filename);
else
    [X,genelist]=i_read_exprmat(filename);
end
%     catch        
%         fprintf('An Expression Matrix file looks like:\n');
%         fprintf('===================\n');
%         fprintf('Gene,c1,c2,c3,c4,c5\n');
%         fprintf('HES4,45,50,15,19,50\n');
%         fprintf('ISG15,168,279,312,425,180\n');
%         fprintf('AGRN,2,3,4,9,5\n');
%         fprintf('C1orf159,0,0,0,1,1\n');
%         fprintf('===================\n');        
%         error('Error hapens when reading Expression Matrix file.');
%     end

end

function [X,genelist,sampleid]=i_read_exprmat(filename,genecolnum,verbose)
% Validate input args
narginchk(1,Inf);

% Get Filename
if ~ischar(filename) && ~(isstring(filename) && isscalar(filename))
    error(message('MATLAB:csvread:FileNameMustBeString')); 
end
filename = char(filename);

% Make sure file exists
if exist(filename,'file') ~= 2 
    error(message('MATLAB:csvread:FileNotFound'));
end


if nargin<2, genecolnum=1; end
if nargin<3, verbose=true; end
if verbose, fprintf('Reading %s ...... ',filename); end
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

T=readtable(filename,'filetype','text');


warning('on','MATLAB:table:ModifiedAndSavedVarnames');
X=table2array(T(:,(1+genecolnum):end));
genelist=string(table2array(T(:,1:genecolnum)));
if nargout>2
    sampleid=string(T.Properties.VariableNames(1+genecolnum:end)');
end
if verbose, fprintf('done.\n'); end
end