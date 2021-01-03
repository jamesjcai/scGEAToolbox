% https://github.com/jblocher/firm_networks/blob/master/matlab/getBonacichCentrality.m
%
% Returns a vector of size n (number of nodes) of argument adjacency matrix
% adj. Input beta is the Bonacich beta factor, and should be in [-1,1].
% Most usages will have beta in [0,1]. If beta = 0, then this is identical
% to degree centrality, if beta = 1 then it is identical to eigenvector
% centrality (to a proportionality constant). Negative betas allow for
% valued relations (positive and negative) - see documentation for more
% info. 
%
% INPUTs: adjacency matrix (required), beta factor (optional. Set to 1 if omitted)
% OUTPUTs: bonacich centrality vector
% Jesse Blocher April 2010
function x=getBonacichCentrality(adj, beta)
if ~isconnected(adj)
    error('Adj matrix not a single component in function bonacich_centrality');
end
if issparse(adj)
    D=eigs(adj);
else
    D=eig(adj);
end
max_eig=max(D);

% Katz-Bonacich Centrality (Beta Centrality)
n = length(adj); % number of nodes
if nargin < 2
    beta = 1; % If 1, should be eigenvector centrality to a constant. 
end;
if beta > 1
    warn('Beta cannot be greater than 1. Setting to 1.');
    beta = 1;
end;
if beta < -1
    warn('Beta cannot be less than -1. Setting to -1.');
    beta = -1;
end;

x = (eye(n)-(beta/max_eig)*adj)\ones(n,1);
