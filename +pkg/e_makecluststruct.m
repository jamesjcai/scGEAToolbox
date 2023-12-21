function [out] = e_makecluststruct
     % see also: e_makeembedstruct
out = struct('kmeans', [], 'snndpc', [], 'sc3', []);
out = orderfields(out);
