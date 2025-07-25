function [edgeInd, gce]=gridEdgeInd(clusterId, M, mins, deltas, pointers)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    cluInd=find(pointers==clusterId);
    gce=edu.stanford.facs.swing.GridClusterEdge(M);
    gce.computeAll(cluInd, mins, deltas)
    edgeInd=gce.edgeBins;
end