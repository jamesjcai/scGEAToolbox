function fname = m2gexf(G, filename)
%M2GEXF Export a MATLAB graph/digraph to GEXF format (Gephi).

if nargin < 2
    fname = [tempname, '.gexf'];
else
    if ~endsWith(filename, '.gexf')
        fname = [filename, '.gexf'];
    else
        fname = filename;
    end
end

if isa(G, 'digraph')
    edgetype = 'directed';
else
    edgetype = 'undirected';
end

ndname = string(G.Nodes.Name);
n = G.numnodes;
A = adjacency(G, 'weighted');
[ii, jj, ww] = find(A);

% Rescale weights to [1, 10] so all edges are visible in Gephi
ww = rescale(abs(ww), 1, 10) .* sign(ww);

fid = fopen(fname, 'w');
fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid, '<gexf xmlns="http://gexf.net/1.2" version="1.2">\n');
fprintf(fid, '  <meta>\n');
fprintf(fid, '    <creator>scGEAToolbox</creator>\n');
fprintf(fid, '  </meta>\n');
fprintf(fid, '  <graph defaultedgetype="%s" mode="static">\n', edgetype);

% Nodes
fprintf(fid, '    <nodes>\n');
for k = 1:n
    lbl = char(ndname(k));
    lbl = strrep(lbl, '&', '&amp;');
    lbl = strrep(lbl, '<', '&lt;');
    lbl = strrep(lbl, '>', '&gt;');
    lbl = strrep(lbl, '"', '&quot;');
    fprintf(fid, '      <node id="%d" label="%s"/>\n', k - 1, lbl);
end
fprintf(fid, '    </nodes>\n');

% Edges
fprintf(fid, '    <edges>\n');
for k = 1:length(ii)
    fprintf(fid, '      <edge id="%d" source="%d" target="%d" weight="%.6g"/>\n', ...
        k - 1, ii(k) - 1, jj(k) - 1, ww(k));
end
fprintf(fid, '    </edges>\n');

fprintf(fid, '  </graph>\n');
fprintf(fid, '</gexf>\n');
fclose(fid);
end
