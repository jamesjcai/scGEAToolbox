% function cl = mcla(cls,k)
%
% DESCRIPTION
%  Performs MCLA for CLUSTER ENSEMBLES
%
% Copyright (c) 1998-2011 by Alexander Strehl

function cl = mcla(cls, k)

disp('CLUSTER ENSEMBLES using MCLA');

if ~exist('k')
    k = max(max(cls));
end

disp('mcla: preparing graph for meta-clustering');
    clb = clstoclbs(cls);
    if false, simbjac(1); end
    cl_lab = clcgraph(clb, k, 'simbjac');
    for i = 1:max(cl_lab)
        matched_clusters = find(cl_lab == i);
        clb_cum(i, :) = mean(clb(matched_clusters, :), 1);
    end
    cl = clbtocl(clb_cum);
