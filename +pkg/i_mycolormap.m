function [colors] = i_mycolormap(n)
if nargin < 1, n = 41; end
plotColor = [; ...
    [.65, .65, .65]; ... % light gray         (0)
[0.1, 0.74, 0.95]; ... % deep sky-blue     (1)
[0.95, 0.88, 0.05]; ... % gold/yellow       (2)
[0.80, 0.05, 0.78]; ... % magenta           (3)
[0.3, 0.8, 0.20]; ... % lime green        (4)
[0.95, 0.1, 0.1]; ... % crimson red       (5)
[0.64, 0.18, 0.93]; ... % blue-violet       (6)
[0.88, 0.56, 0]; ... % orange            (7)
[0.4, 1.0, 0.7]; ... % aquamarine        (8)
[0.95, 0.88, 0.7]; ... % salmon-yellow     (9)
[0, 0.2, 1]; ... % blue              (10)
[1, 0.41, 0.7]; ... % hot pink          (11)
[0.5, 1, 0]; ... % chartreuse        (12)
[0.6, 0.39, 0.8]; ... % amtheyist         (13)
[0.82, 0.36, 0.36,]; ... % indian red        (14)
[0.53, 0.8, 0.98]; ... % light sky blue    (15)
[0, 0.6, 0.1]; ... % forest green      (16)
[0.65, 0.95, 0.5]; ... % light green       (17)
[0.85, 0.6, 0.88]; ... % light purple      (18)
[0.90, 0.7, 0.7]; ... % light red         (19)
[0.2, 0.2, 0.6]; ... % dark blue         (20)
];

repeats = max(1, ceil(n/19));
colors = [plotColor(1, :); repmat(plotColor(2:end, :), repeats, 1)];
end
