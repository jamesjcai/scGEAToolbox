function [x_centre, y_centre, r, sq_error, th] = circleFit ( P )
%# CIRCFIT fits a circle to a set of points using least sqaures
%#  P is a 2 x n matrix of points to be fitted

per_error = 0.1/100; % i.e. 0.1%

%# initial estimates
X  = mean(P, 2)';
r = sqrt(mean(sum((repmat(X', [1, length(P)]) - P).^2)));

v_cen2points = zeros(size(P));
niter = 0;

%# looping until convergence
while niter < 1 || per_diff > per_error

    %# vector from centre to each point
    v_cen2points(1, :) = P(1, :) - X(1);
    v_cen2points(2, :) = P(2, :) - X(2);  

    %# distacnes from centre to each point
    centre2points = sqrt(sum(v_cen2points.^2));

    %# distances from edge of circle to each point
    d = centre2points - r;

    %# computing 3x3 jacobean matrix J, and solvign matrix eqn.
    R = (v_cen2points ./ [centre2points; centre2points])';
    J = [ -ones(length(R), 1), -R ];
    D_rXY = -J\d';

    %# updating centre and radius
    r_old = r;    X_old = X;
    r = r + D_rXY(1);
    X = X + D_rXY(2:3)';

    %# calculating maximum percentage change in values
    per_diff = max(abs( [(r_old - r) / r, (X_old - X) ./ X ])) * 100;

    %# prevent endless looping
    niter = niter + 1;
    if niter > 1000
        %break;
        error('Convergence not met in 1000 iterations!')
    end
end

x_centre = X(1);
y_centre = X(2);
sq_error = sum(d.^2);

end


