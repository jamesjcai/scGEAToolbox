function e_m2cyto(obj, filename_in)
% M2CYTO - Converting @graph object to Cytoscape-readable format.
% function [filename] = m2cyto(obj,filename_in,link_attrib)

% WRITTEN BY       : Kenth Engø-Monsen, 2008.12.12
% LAST MODIFIED BY : Kenth Engø-Monsen, 2012.04.18

A = obj.adjacency; % Adjacency matrix for the graph object.
n     = obj.numnodes;

link_attrib = obj.Edges.Weight;
ndname = string(obj.Nodes.Name);

if nargin > 1
    fname = [filename_in, '.m2c'];
else
    fname = [inputname(1), '.m2c'];
end

fid = fopen(fname, 'w');
[ii, jj, ~] = find(A);
[~, ~, att] = find(link_attrib);
fprintf(fid, 'Source\tTarget\tWeight\n');
for ik = 1:length(att) % nnz(A)
    fprintf(fid, '%s\t%s\t%f\n', ndname(ii(ik)), ndname(jj(ik)), att(ik));
end
fclose(fid);
end
