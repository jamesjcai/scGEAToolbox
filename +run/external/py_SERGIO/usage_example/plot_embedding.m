function plot_embedding(x, coords, coloring, labels)

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
% x - adjacency matrix (NxN) of the network
%
% coords - polar (Nx2) or spherical (Nx3) hyperbolic coordinates of the nodes
%   in the hyperbolic disk they are in the form: [theta,r]
%   in the hyperbolic sphere they are in the form: [azimuth,elevation,r]
%
% coloring - string indicating how to color the nodes:
%   'popularity' - nodes colored by degree with a blue-to-red colormap
%                  (valid for 2D and 3D)
%   'similarity' - nodes colored by angular coordinate with a HSV colormap
%                  (valid only for 2D)
%   'labels' - nodes colored by labels, which can be all unique
%              (for example to indicate an ordering of the nodes)
%              or not (for example to indicate community memberships)
%              (valid for 2D and 3D)
%
% labels - numerical labels for the nodes (only needed if coloring = 'labels')

% check input
narginchk(3,4);
validateattributes(x, {'numeric'}, {'square','finite','nonnegative'});
if ~issymmetric(x)
    error('The input matrix must be symmetric.')
end
if any(x(speye(size(x))==1))
    error('The input matrix must be zero-diagonal.')
end
validateattributes(coords, {'numeric'}, {'2d','nrows',length(x)})
dims = size(coords,2);
validateattributes(dims, {'numeric'}, {'>=',2,'<=',3});
validateattributes(coloring, {'char'}, {});
if dims == 2 && ~any(strcmp(coloring,{'popularity','similarity','labels'}))
    error('Possible coloring options in 2D: ''popularity'',''similarity'',''labels''.');
end
if dims == 3 && ~any(strcmp(coloring,{'popularity','labels'}))
    error('Possible coloring options in 3D: ''popularity'',''labels''.');
end
if strcmp(coloring,'labels')
    validateattributes(labels, {'numeric'}, {'vector','numel',length(x)})
end

% set plot options
edge_width = 1;
edge_color = [0.85 0.85 0.85];
node_size = 150;

% set the node colors
if strcmp(coloring,'popularity')
    deg = full(sum(x>0,1));
    deg = round((max(deg)-1) * (deg-min(deg))/(max(deg)-min(deg)) + 1);
    colors = colormap_blue_to_red(max(deg));
    colors = colors(deg,:);
elseif strcmp(coloring,'similarity')
    colormap('hsv')
    colors = coords(:,1);
elseif strcmp(coloring,'labels')
    uniq_lab = unique(labels);
    temp = zeros(size(labels));
    for i = 1:length(uniq_lab)
        temp(labels==uniq_lab(i)) = i;
    end
    labels = temp; clear uniq_lab temp;
    colors = hsv(length(unique(labels)));
    colors = colors(labels,:);
end

% plot the network
hold on
radius = 2*log(length(x));
if dims == 2
    [coords(:,1),coords(:,2)] = pol2cart(coords(:,1),coords(:,2));
    [h1,h2] = gplot(x, coords, 'k'); plot(h1, h2, 'Color', edge_color, 'LineWidth', edge_width);
    scatter(coords(:,1), coords(:,2), node_size, colors, 'filled', 'MarkerEdgeColor', 'k');
    xlim([-radius, radius]); ylim([-radius, radius])
elseif dims == 3
    [coords(:,1),coords(:,2),coords(:,3)] = sph2cart(coords(:,1),coords(:,2),coords(:,3));
    [r,c] = find(triu(x>0,1));
    for i = 1:length(r)
        plot3([coords(r(i),1) coords(c(i),1)],[coords(r(i),2) coords(c(i),2)], [coords(r(i),3) coords(c(i),3)], ...
            'Color', edge_color, 'LineWidth', edge_width)
    end
    scatter3(coords(:,1), coords(:,2), coords(:,3), node_size, colors, 'filled', 'MarkerEdgeColor', 'k');
    xlim([-radius, radius]); ylim([-radius, radius]); zlim([-radius, radius])
end
axis square
axis off

function colors = colormap_blue_to_red(n)

colors = zeros(n,3);
m = round(linspace(1,n,4));
colors(1:m(2),2) = linspace(0,1,m(2));
colors(1:m(2),3) = 1;
colors(m(2):m(3),1) = linspace(0,1,m(3)-m(2)+1);
colors(m(2):m(3),2) = 1;
colors(m(2):m(3),3) = linspace(1,0,m(3)-m(2)+1);
colors(m(3):n,1) = 1;
colors(m(3):n,2) = linspace(1,0,n-m(3)+1);
