function T = e_groupsummary(y, c)
%E_GROUPSUMMARY Per-group mean, SEM and median of one value per cell.
%   T = pkg.e_groupsummary(Y, C) returns a table with one row per group of
%   C and the variables Group, Mean, SEM and Median of Y within it. Rows
%   are in FINDGROUPS order (sorted labels), the order every caller labels
%   its axes and tables with.
%
%   This exists because GRPSTATS orders text groups by first appearance,
%   not sorted: bars and exported values computed with GRPSTATS but
%   labelled from FINDGROUPS were attached to the wrong groups whenever the
%   cells did not already come in sorted order.
%
%   Y - numeric vector, one value per cell
%   C - grouping vector of the same length (string, cellstr, categorical
%       or numeric)
%
%   See also findgroups, splitapply.

arguments
    y {mustBeNumeric, mustBeVector}
    c {mustBeVector}
end

if numel(y) ~= numel(c)
    error("pkg:e_groupsummary:SizeMismatch", ...
        "Y has %d values but C has %d; they must be one per cell.", ...
        numel(y), numel(c));
end

y = double(y(:));
[gid, Group] = findgroups(string(c(:)));
Mean = splitapply(@mean, y, gid);
SEM = splitapply(@(v) std(v)/sqrt(numel(v)), y, gid);
Median = splitapply(@median, y, gid);
T = table(Group, Mean, SEM, Median);
end
