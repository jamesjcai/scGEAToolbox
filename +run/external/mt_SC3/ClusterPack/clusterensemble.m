% function cl = clusterensemble(cls,k)
%
% DESCRIPTION
%   - performs the supra-consensus function for CLUSTER ENSEMBLES
%     to combine multiple clustering label vectors e.g., cls = [cl1; cl2;...]
%   - returns a row vector of the combined clustering labels
%   - the following consensus functions are computed
%     - Cluster-based Similiarty Partitioning Algorithm (CSPA)
%     - HyperGraph Partitioning Algorithm (HGPA)
%     - Meta-CLustering Algorithm (MCLA)
%     and the one with the maximum average normalized mutual information
%     is returned
% ARGUMENTS
%   cls   - matrix of labelings, one labeling per row (n x p)
%           entries must be integers from 1 to k1 (row 1), 1 to k2 (row 2),...
%           use NaN as an entry for unknown/missing labels
%   k     - 1,2,3,... maximum number of clusters in the combined clustering
%           (optional, default max(max(cls))
% EXAMPLE
%   clusterensemble;
%   clusterensemble([ones(3,20) 2*ones(3,30); ones(1,10) 2*ones(1,40)])
% REFERENCE
%   please refer to the following paper if you use CLUSTER ENSEMBLES
%     A. Strehl and J. Ghosh. "Cluster Ensembles - A Knowledge Reuse
%	  Framework for Combining Multiple Partitions", Journal on
%	  Machine Learning Research (JMLR), 3:583-617, December 2002
% RELEASE
%   version 2.0, 2011/05/01, tested on WIN7 Octave 3.2.4 and LNX86 Matlab 5.2.0.3084
%   available from http://www.strehl.com
%   license granted for research use ONLY (see README)
%   copyright (c) 1998-2011 by Alexander Strehl

function cl = clusterensemble(cls, k)

if false
    cls = [1, 1, 1, 2, 2, 3, 3; 2, 2, 2, 3, 3, 1, 1; 1, 1, 2, 2, 3, 3, 3; 1, 2, NaN, 1, 2, NaN, NaN];
    cl = mcla(cls, 2);
    cl = hgpa(cls, 2);
    cl = cspa(cls, 2);
end

if ~exist('cls', 'var')
    disp('clusterensemble-warning: no arguments - displaying illustrative example:');
        cls = [1, 1, 1, 2, 2, 3, 3; 2, 2, 2, 3, 3, 1, 1; 1, 1, 2, 2, 3, 3, 3; 1, 2, NaN, 1, 2, NaN, NaN];
        disp('clusterensemble-advice: type "help clusterensemble" for information about usage');
            disp(' ');
        end

        if size(cls, 2) > 1000
            workfcts = {'hgpa', 'mcla'};
            disp('clusterensemble-warning: using only hgpa and mcla');
        else
            workfcts = {'cspa', 'hgpa', 'mcla'};
        end

        for i = 1:length(workfcts)
            workfct = workfcts{i};
            if ~exist('k', 'var')
                cl(i, :) = feval(workfct, cls);
            else
                cl(i, :) = feval(workfct, cls, k);
            end
            q(i) = ceevalmutual(cls, cl(i, :));
            disp(['clusterensemble: ', workfct, ' at ', num2str(q(i))]);
        end

        [~, best] = max(q);

        cl = cl(best, :);
