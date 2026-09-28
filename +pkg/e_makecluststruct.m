function [out] = e_makecluststruct
     % see also: e_makeembedstruct
out = struct('kmeans', [], 'snndpc', [], 'louvain', [], 'louvainpc', [], 'sc3', []);
out = orderfields(out);
