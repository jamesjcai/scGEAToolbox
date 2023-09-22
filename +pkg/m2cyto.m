function m2cyto(obj, filename_in)
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
%if nargin == 3,
[~, ~, att] = find(link_attrib);
fprintf(fid, 'Source\tTarget\tWeight\n');
for ik = 1:length(att) % nnz(A)
    fprintf(fid, '%s\t%s\t%f\n', ndname(ii(ik)), ndname(jj(ik)), att(ik));
    %     fprintf(fid,[num2str(ii(i)) ' ' num2str(jj(i)) ' ' num2str(val(i)) ' ' ...
    %                  num2str(att(i)) '\n']);
end
% else
%   for i = 1:nnz(m)
%     fprintf(fid,[num2str(ii(i)) ' ' num2str(jj(i)) ' ' num2str(val(i)) '\n']);
%   end
%   fprintf(fid,[num2str(n) ' ' num2str(n) ' ' num2str(0) '\n']);
% end
fclose(fid);
end
