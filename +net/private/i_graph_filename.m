function fname = i_graph_filename(filename, ext)
% I_GRAPH_FILENAME  Resolve output filename for graph export functions.
%   fname = i_graph_filename(filename, ext) returns filename with the
%   correct extension. If filename is empty, a temporary file is used.
%
%   Example:
%     fname = i_graph_filename('mynet', 'graphml')   % -> 'mynet.graphml'
%     fname = i_graph_filename('', 'cx')             % -> '<tempname>.cx'
if isempty(filename)
    fname = [tempname, '.', ext];
elseif ~endsWith(filename, ['.', ext])
    fname = [filename, '.', ext];
else
    fname = filename;
end
end
