function [coords_2_aligned, aligned_corr] = angular_alignment(coords_1, ...
    coords_2, rotations)

% Authors:
% - main code: Alessandro Muscoloni, 2017-09-21

% Released under MIT License
% Copyright (c) 2017 A. Muscoloni, J. M. Thomas, C. V. Cannistraci

% Reference:
% A. Muscoloni, J. M. Thomas, S. Ciucci, G. Bianconi, and C. V. Cannistraci,
% "Machine learning meets complex networks via coalescent embedding in the hyperbolic space",
% Nature Communications 8, 1615 (2017). doi:10.1038/s41467-017-01825-5

%%% INPUT %%%
% coords_1  - angular coordinates 1
% coords_2  - angular coordinates 2
% rotations - number of angular rotations to apply

%%% OUTPUT %%%
% coords_2_aligned - angular coordinates 2 aligned to angular coordinates 1
%                    maximizing the correlation between them
% aligned_corr     - correlation maximized by the alignment

% force to use column vectors
if ~iscolumn(coords_1)
    coords_1 = coords_1';
end
if ~iscolumn(coords_2)
    coords_2 = coords_2';
end

% initialization
coords_2_shifted = cell(2*rotations,1);
shifted_corr = zeros(2*rotations,1);

% reverse the angular coordinates 2
coords_2_rev = abs(coords_2 - 2*pi);

for i = 1:rotations
    % apply "i" rotations to the angular coordinates 2
    coords_2_shifted{i} = mod(coords_2 + i*(2*pi/rotations), 2*pi);
    % apply "i" rotations to the reversed angular coordinates 2
    coords_2_shifted{rotations+i} = mod(coords_2_rev + i*(2*pi/rotations), 2*pi);
    % compute the correlations
    shifted_corr(i) = corr(coords_1, coords_2_shifted{i});
    shifted_corr(rotations+i) = corr(coords_1, coords_2_shifted{rotations+i});
end

% retrieve the shifted angular coordinates 2
% that maximize the correlation
[aligned_corr, idx] = max(shifted_corr);
coords_2_aligned = coords_2_shifted{idx};