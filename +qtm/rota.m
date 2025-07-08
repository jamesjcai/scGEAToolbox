function rotated_coords = rota(x, y, z, alpha, beta, gamma)
% ROTATE_POINT Rotates a 3D point (x, y, z) using extrinsic Z-Y-X rotation
% alpha, beta, gamma are in radians

% Full rotation matrix (extrinsic Z -> Y -> X)
R = [ cos(alpha)*cos(beta),  cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma),  cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma);
      sin(alpha)*cos(beta),  sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma),  sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma);
     -sin(beta),             cos(beta)*sin(gamma),                                     cos(beta)*cos(gamma) ];

% Original point vector
v = [x; y; z];

% Apply rotation
rotated_coords = R * v;

% Display result
fprintf('Original point:  (%.4f, %.4f, %.4f)\n', x, y, z);
fprintf('Rotated point:   (%.4f, %.4f, %.4f)\n', rotated_coords(1), rotated_coords(2), rotated_coords(3));
end