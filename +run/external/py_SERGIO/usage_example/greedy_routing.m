function [GR_score, succ_paths, avg_length] = greedy_routing(x, coords)

% Authors:
% - main code: Alessandro Muscoloni, 2017-09-21

% Released under MIT License
% Copyright (c) 2017 A. Muscoloni, J. M. Thomas, C. V. Cannistraci

% Reference:
% A. Muscoloni, J. M. Thomas, S. Ciucci, G. Bianconi, and C. V. Cannistraci,
% "Machine learning meets complex networks via coalescent embedding in the hyperbolic space",
% Nature Communications 8, 1615 (2017). doi:10.1038/s41467-017-01825-5

% Note: the function contains parallel loops to speedup the computation.

% Greedy routing:
% at each step the packet is sent to the node neighbour closest to the
% destination in the hyperbolic space.
% the packet is dropped when a loop is detected.

%%% INPUT %%%
% x - adjacency matrix (NxN) of the network
%
% coords - polar (Nx2) or spherical (Nx3) hyperbolic coordinates of the nodes
%   in the hyperbolic disk they are in the form: [theta,r]
%   in the hyperbolic sphere they are in the form: [azimuth,elevation,r]

%%% OUTPUT %%%
% GR_score   - score between 0 and 1
%              average ratio between shortest paths and hop-length of GR paths
%              the hop-length of unsuccessful GR paths is infinite
%
% succ_paths - percentage of successful GR paths
%
% avg_length - average hop-length of the successful GR paths

% check input
narginchk(2,2);
validateattributes(x, {'numeric'}, {'square','finite','nonnegative'});
if ~issymmetric(x)
    error('The input matrix must be symmetric.')
end
if any(x(speye(size(x))==1))
    error('The input matrix must be zero-diagonal.')
end
n = size(x,1);
validateattributes(coords, {'numeric'}, {'2d','nrows',n})
dims = size(coords,2);
validateattributes(dims, {'numeric'}, {'>=',2,'<=',3});

% hyperbolic pairwise distances
if dims == 2
    hdist = squareform(pdist(coords, @hyperbolic_dist_2D));
elseif dims == 3
    hdist = squareform(pdist(coords, @hyperbolic_dist_3D));
end

% for each pair (i,j), compute the neighbour of i closest to j
next = zeros(size(x));
parfor j = 1:n
    next(:,j) = compute_next_to_j(j, x, hdist(:,j), n);
end

% paths contains:
% Inf    -> unsuccessful path
% >= 0   -> hop-length of successful path
paths = zeros(size(x));

% for each pair (i,j), greedy routing from node i to node j
parfor j = 1:n
    paths(:,j) = greedy_routing_to_j(j, next(:,j), n);
end

% percentage of successful paths (note: the diagonal is not considered)
succ_paths = length(paths(paths~=0 & paths~=Inf)) / length(paths(paths~=0));

% average hop-length of the successful paths
avg_length = mean(paths(paths~=0 & paths~=Inf));

% GR-score
%sh_paths = graphallshortestpaths(sparse(x), 'Directed', false);
sh_paths = distances(graph(sparse(x)));

GR_score = mean(sh_paths(paths~=0) ./ paths(paths~=0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function next_to_j = compute_next_to_j(j, x, hdist_to_j, n)

next_to_j = zeros(n,1);
for i = 1:n
    if i ~= j
        temp = hdist_to_j;
        temp(x(i,:)==0) = Inf;
        [~,next_to_j(i)] = min(temp);
    end
end

function paths_to_j = greedy_routing_to_j(j, next_to_j, n)

paths_to_j = zeros(n,1);
for i = 1:n
    if i ~= j && paths_to_j(i) == 0
        paths_to_j = greedy_routing_rec(i, j, next_to_j, paths_to_j);
    end
end

function paths_to_j = greedy_routing_rec(i, j, next_to_j, paths_to_j)

if next_to_j(i) == j
    paths_to_j(i) = 1;
elseif next_to_j(next_to_j(i)) == i
    paths_to_j(i) = Inf;
else
    if paths_to_j(next_to_j(i)) == 0
        paths_to_j = greedy_routing_rec(next_to_j(i), j, next_to_j, paths_to_j);
    end
    if paths_to_j(next_to_j(i)) == Inf
        paths_to_j(i) = Inf;
    else
        paths_to_j(i) = 1 + paths_to_j(next_to_j(i));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = hyperbolic_dist_2D(XI, XJ)
%Computes the hyperbolic distance between point XI (1 x 2) and the m points
%stored in XJ (m x 2). The coordinates of these points are polar in the
%format (angular coord, radial coord). The resulting similarities are
%stored in d.
%
%INPUT
%   XI -> The polar coordinates of a single point in the Poincaré disc.
%   XJ -> The polar coordinates of m points in the Poincaré disc.
%
%OUTPUT
%   d -> The hyperbolic distance between point XI and the other m points
%        stored in XJ. The hyperbolic distance between points (Ti, Ri) and
%        (Tj, Rj) in the hyperbolic space H^2 of curvature K = -1,
%        represented by the Poincaré disc is:
%
% Dij = arccosh(cosh(Ri)*cosh(Rj) - sinh(Ri)*sinh(Rj)*cos(Tij));
%
%        with Tij = pi - |pi - |Ti - Tj||
%
% Copyright (C) Gregorio Alanis-Lobato, 2014

A =  pi - abs(pi - abs(XI(1) - XJ(:,1))); %angular separation between points
d = acosh(cosh(XI(2)).*cosh(XJ(:,2)) - sinh(XI(2)).*sinh(XJ(:,2)).*cos(A));
d(isinf(d)) = 0;

% due to numerical approximations, points with a tiny or zero angular
% separation and close radial coordinates could produce a wrong complex
% hyperbolic distance with zero real part. These distances are replaced
% by the radial separation as expected by the theoretical formula.
if ~isreal(d)
   d(imag(d)~=0) = abs(XI(2)-XJ(imag(d)~=0,2)); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = hyperbolic_dist_3D(XI, XJ)

% Hyperbolic distance between points in the Poincaré model of 3D hyperbolic space
% with curvature K = -1.
% The function can be used as custom distance for the pdist function.
% The coordinates are in spherical format (azimuth, elevation, radius).

%%% INPUT %%%
% XI - (1 x 3) coordinates of one point
% XJ - (m x 3) coordinates of m points

%%% OUTPUT %%%
% D - (m x 1) distances between the point XI and the m points XJ

% compute the angular separation between the points
% note: it is the angle between the two radii going from the center of the
% sphere to the points; since the radii are straight lines, this angle
% is the same both in the euclidean and hyperbolic space; therefore it can
% be obtained computing the euclidean distance between the two points
% and solving a SSS triangle.
CI = zeros(size(XI));
CJ = zeros(size(XJ));
[CI(1),CI(2),CI(3)] = sph2cart(XI(1),XI(2),XI(3));
[CJ(:,1),CJ(:,2),CJ(:,3)] = sph2cart(XJ(:,1),XJ(:,2),XJ(:,3));
A = acos((XI(3)^2 + XJ(:,3).^2 - sum((repmat(CI,size(CJ,1),1)-CJ).^2,2)) ./ (2*XI(3)*XJ(:,3)));
A(isnan(A) | isinf(A) | ~isreal(A)) = 0; % degenerate case (points in the center)

% compute the hyperbolic distance applying the hyperbolic law of cosines
D = acosh(cosh(XI(3)).*cosh(XJ(:,3)) - sinh(XI(3)).*sinh(XJ(:,3)).*cos(A));

% due to numerical approximations, points with a tiny or zero angular
% separation and close radial coordinates could produce a wrong complex
% hyperbolic distance with zero real part. These distances are replaced
% by the radial separation as expected by the theoretical formula.
if ~isreal(D)
   D(imag(D)~=0) = abs(XI(3)-XJ(imag(D)~=0,3)); 
end