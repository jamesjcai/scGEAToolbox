function [A, s, G] = genie3(X, genelist, donorm, plotit)
% GENIE3: GEne Network Inference with Ensemble of trees
%
% [A] = net.genie3(X)
% [A, s, G] = net.genie3(X, genelist, donorm, plotit)
%
% Uses random forest regression to infer gene regulatory network.
% ref: Huynh-Thu et al., PLoS ONE, 2010.

if nargin < 2 || isempty(genelist)
    genelist = string(1:size(X, 1));
end
if nargin < 3, donorm = false; end
if nargin < 4, plotit = false; end

pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '..', 'external', 'ml_GENIE3');
if ~(ismcc || isdeployed), addpath(pth); end
pth = fullfile(pw1, '..', 'external', 'ml_GENIE3', 'RT');
if ~(ismcc || isdeployed), addpath(pth); end

if donorm
    X = sc_norm(X, 'type', 'libsize');
end
% GENIE3 takes data (rows: cell, columns: genes)
data = full(X.');

A = GENIE3_ori(data);
if nargout == 2
    s = get_link_list(A);
end

%%
if plotit || nargout == 3
    s = get_link_list(A);
    G = digraph(s(:, 1), s(:, 2), s(:, 3), genelist);
end
if plotit
    LWidths = 5 * G.Edges.Weight / max(G.Edges.Weight);
    LWidths(LWidths == 0) = 1e-5;
    try
        p = plot(G, 'LineWidth', LWidths);
        p.MarkerSize = 7;
        p.Marker = 's';
        p.NodeColor = 'r';
    catch ME
        warning(ME.message);
    end
end
end
