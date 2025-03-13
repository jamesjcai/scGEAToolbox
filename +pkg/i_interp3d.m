function [xyz2] = i_interp3d(xyz, numInterpPoints)

if nargin<2, numInterpPoints = size(xyz,1)*2; end

distances = sqrt(sum(diff(xyz).^2, 2)); % Euclidean distances between points
arcLength = [0; cumsum(distances)]; % Arc length parameter

% Define new parameter values for interpolation
% numInterpPoints = 1000;
arcLengthInterp = linspace(0, arcLength(end), numInterpPoints); 

% Interpolate x, y, and z separately
xInterp = interp1(arcLength, xyz(:,1), arcLengthInterp, 'spline');
yInterp = interp1(arcLength, xyz(:,2), arcLengthInterp, 'spline');
zInterp = interp1(arcLength, xyz(:,3), arcLengthInterp, 'spline');

xyz2 = [xInterp' yInterp' zInterp'];
end

%{
% Generate example 3D curve (replace with your actual data)
numPoints = 46;
t = linspace(0, 2*pi, numPoints)'; % Example parametric variable
x = cos(t);  % Example x-coordinates
y = sin(t);  % Example y-coordinates
z = t;       % Example z-coordinates

% Compute cumulative arc length
xyz = [x, y, z]; 

% Plot results
figure;
plot3(x, y, z, 'ro', 'MarkerFaceColor', 'r'); hold on; % Original points
plot3(xInterp, yInterp, zInterp, 'b-', 'LineWidth', 1.5); % Interpolated curve
grid on; xlabel('X'); ylabel('Y'); zlabel('Z');
legend('Original 46 points', 'Interpolated 1000 points');
title('3D Curve Interpolation');
axis equal;
%}