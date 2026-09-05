function [x, y] = quadraticcurveto(currentp, curvetop, controllp)
%QUADRATICCURVETO Points along a quadratic Bezier between two points.
%
%   Lifted from the local copy in gui.i_networkvis so that
%   gui.i_networkvisccc can draw the same curves without a third copy of it.
%   gui.i_networkvis and gui.i_networkvis_linear still carry their own
%   locals, which take precedence over this file, so adding it here changes
%   nothing for them.
%
%   Same maths, two tidies: the output is preallocated, and three variables
%   the original assigned but never read (x1, y1, x2 - used only by an
%   alternative control point that is commented out) are gone.

    % The quadraticCurveTo method creates a line from the path's current point to the specified point, via a controlpoint.
    % The quadraticCurveTo method takes four parameters, the x and y coordinates for the controlpoint, and the x and y coordinates for the line's destination.

    % Current point moveTo(20,20)
    % CurveTo point quadraticCurveTo(20,100,200,20)
    % Controllpoint quadraticCurveTo(20,100,200,20)

    o = curvetop(1);
    i = curvetop(2);
    d = currentp(1);
    f = currentp(2);

    if nargin < 3
        controllp = [(d + o) / 2 + (i - f) / 4, (f + i) / 2 + (d - o) / 4];
    end
    P = [[currentp'; 0], [controllp'; 0], [curvetop'; 0]];

    n = 3;
    count = 1;
    div = 50; % number of segments of the curve (Increase this value to obtain a
    % smoother curve
    A = zeros(3, div + 1);
    for u = 0:(1 / div):1
        sum = [0, 0, 0]';
        for i = 1:n
            B = nchoosek(n, i-1) * (u^(i - 1)) * ((1 - u)^(n - i + 1)); % B is the Bernstein polynomial value
            sum = sum + B * P(:, i);
        end
        B = nchoosek(n, n) * (u^(n));
        sum = sum + B * P(:, n);
        A(:, count) = sum; % the matrix containing the points of curve as column vectors.
        count = count + 1; % count is the index of the points on the curve.
    end

    x = A(1, :);
    y = A(2, :);
 end
