function [eigenvalues,No_cluster,No_cluster2] = Num_cluster(W,NC1)
% Estimation of the number of clusters from single cell data
%
% Input:
%   -- W: cell-to-cell similarity matrix
%   -- NC1: number of clusters specified
%
% Output:
%   -- eigenvalues: eigenvalues of the graph Laplacian of the truncated
%   consensus matrix.
%   -- No_cluster: Number of clusters inferred by SoptSC.
%
nno_cluster = 2:20;  % 2:20

if NC1 <= 5
    tau = 0.3; % 0.4
elseif NC1 <= 10
    tau = 0.4;
else
    tau = 0.5;
end

tol = 0.01;
[all_eigs,~] = consen(W,nno_cluster,tau);

zz = sort(abs(real(all_eigs)));

if length(zz)>=2
    gap = zz(2:end) - zz(1:end-1);
    [~,No_cluster1] = max(gap);
end

No_cluster2 = length(find(zz<=tol));
No_cluster = No_cluster1;
display('Number of cluster based on zero eigenvalues & Largest gap ');
display([No_cluster2 No_cluster]);

eigenvalues = zz;
end