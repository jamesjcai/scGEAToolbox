function [HD_corr, C_score] = compare_embedding(coords_1, coords_2)

% Authors:
% - main code: Alessandro Muscoloni, 2017-09-21
% - support functions: indicated at the beginning of the function

% Released under MIT License
% Copyright (c) 2017 A. Muscoloni, J. M. Thomas, C. V. Cannistraci

% Reference:
% A. Muscoloni, J. M. Thomas, S. Ciucci, G. Bianconi, and C. V. Cannistraci,
% "Machine learning meets complex networks via coalescent embedding in the hyperbolic space",
% Nature Communications 8, 1615 (2017). doi:10.1038/s41467-017-01825-5

%%% INPUT %%%
% coords_1, coords_2 - two hyperbolic embeddings to compare:
%   polar (Nx2) or spherical (Nx3) hyperbolic coordinates of the nodes
%   in the hyperbolic disk they are in the form: [theta,r]
%   in the hyperbolic sphere they are in the form: [azimuth,elevation,r]

%%% OUTPUT %%%
% HD_corr - correlation between all pairwise hyperbolic distances
%           (valid for 2D and 3D)
%
% C_score - score evaluating the match between the circular order
%           of the nodes (computed only for 2D)

% check input
narginchk(2,2);
validateattributes(coords_1, {'numeric'}, {'2d'})
validateattributes(coords_2, {'numeric'}, {'2d'})
if size(coords_1,1)~=size(coords_2,1) || size(coords_1,2)~=size(coords_2,2)
    error('Dimensions of arrays coords_1 and coords_2 do not match.')
end
dims = size(coords_1,2);
validateattributes(dims, {'numeric'}, {'>=',2,'<=',3});

% HD-corr
if dims == 2
    HD_1 = pdist(coords_1, @hyperbolic_dist_2D);
    HD_2 = pdist(coords_2, @hyperbolic_dist_2D);
elseif dims == 3
    HD_1 = pdist(coords_1, @hyperbolic_dist_3D);
    HD_2 = pdist(coords_2, @hyperbolic_dist_3D);
end
HD_corr = corr(HD_1', HD_2');

% C-score
C_score = NaN;
if dims == 2
    C_score = compute_c_score(coords_1(:,1), coords_2(:,1));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c_score = compute_c_score(theta_1, theta_2)

%%% INPUT %%%
% theta_1, theta_2 - angular coordinates of two hyperbolic embeddings

%%% OUTPUT %%%
% c_score - percentage of pairs in the same circular order

% note: a pair (i,j) is ordered clockwise if the minimum angular distance
%       from i to j is computed going clockwise (viceversa for counterclockwise).
%       when the distance is 0 or pi, it is treated as a special case

% since the network could be embedded in the opposite clock direction
% both the directions are considered and the maximum is taken
count_dir1 = 0; % correct inference considering one direction of embedding
count_dir2 = 0; % correct inference considering the opposite direction of embedding

n = length(theta_1);
for i = 1:n
    for j = (i+1):n
        
        % shift the pairs so that node i correspond to the origin
        % and then look at the position of node j
        th1_j = mod(theta_1(j) + (2*pi - theta_1(i)), 2*pi);
        th2_j = mod(theta_2(j) + (2*pi - theta_2(i)), 2*pi);
        
        % check pairs arrangement
        if ...      % special cases
                (th1_j==0  && th2_j==0)  || ...
                (th1_j==pi && th2_j==pi)
            count_dir1 = count_dir1 + 1;
            count_dir2 = count_dir2 + 1;
        elseif ...  % networks has same direction
                (th1_j<pi  && th2_j<pi)  || ...
                (th1_j>pi  && th2_j>pi)
            count_dir1 = count_dir1 + 1;
        elseif ...  % networks has opposite direction
                (th1_j<pi  && th2_j>pi)  || ...
                (th1_j>pi  && th2_j<pi)
            count_dir2 = count_dir2 + 1;
        end
    end
end

tot = nchoosek(n,2);    % combinations of pairs
c_score = max(count_dir1/tot, count_dir2/tot);