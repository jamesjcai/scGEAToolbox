function tf = i_isgroupingvar(v, nCells, allowsingle)
%I_ISGROUPINGVAR  True when a per-cell vector can be used to group cells.
%
%   tf = pkg.i_isgroupingvar(v, nCells)
%   tf = pkg.i_isgroupingvar(v, nCells, allowsingle)
%
% Decides whether a named cell attribute is a categorical label, which the
% grouping dialogs can offer, or a per-cell measurement, which they cannot.
% SCE.LIST_CELL_ATTRIBUTES holds both: R_READSEURATRDS stashes meta.data
% columns like seurat_clusters there, and CALLBACK_DRAWTRAJECTORY writes
% pseudotime into the same list. Offering the latter as a grouping variable
% would ask for one group, and one label, per cell.
%
% A variable qualifies when all of these hold:
%   * it has one value per cell
%   * it is a string, cellstr, categorical, logical or numeric vector
%   * if numeric, every value is finite and a whole number -- pseudotime,
%     percent.mt and other measurements are excluded on this test, while
%     integer codes such as seurat_clusters are kept
%   * it has at most MAXLEVELS distinct values, and fewer than one per cell
%     (a value unique to every cell is an identifier, like a barcode)
%
% ALLOWSINGLE false additionally requires more than one level, matching the
% meaning the flag has in GUI.I_SELECT1CLASS: list a variable only when it
% could actually separate the cells.
%
% See also GUI.I_SELECT1CLASS, GUI.I_SELECTNCLASS

arguments
    v
    nCells (1,1) double {mustBeNonnegative}
    allowsingle (1,1) logical = true
end

% Above this, a "grouping" is not one any plot or test could use, and is
% far more likely to be a continuous variable that happens to be integral.
MAXLEVELS = 100;

tf = false;

% A char row vector has one element per character, not per cell, and a char
% matrix is not a per-cell vector either. Cellstr is the text form that is.
if ischar(v), return; end

if ~(isstring(v) || iscellstr(v) || iscategorical(v) || islogical(v) || isnumeric(v))
    return;
end

if isempty(v) || numel(v) ~= nCells, return; end
v = v(:);

if isnumeric(v)
    if ~all(isfinite(v)), return; end
    if any(mod(double(v), 1) ~= 0), return; end
end

nlev = i_countlevels(v);
if nlev < 1, return; end            % nothing but missing values
if nlev > MAXLEVELS, return; end
if nlev == nCells, return; end
if ~allowsingle && nlev <= 1, return; end

tf = true;
end


function n = i_countlevels(v)
%I_COUNTLEVELS Distinct values, not counting the missing ones.
%   UNIQUE treats <undefined> and <missing> the way it treats NaN: each one
%   is its own element. Counting levels with it therefore makes a two-level
%   annotation with 200 unlabeled cells look like 202 levels, which pushes a
%   perfectly good grouping past MAXLEVELS and drops it from the list. An
%   h5ad file that leaves some cells unlabeled is the ordinary case, not an
%   odd one.
if iscategorical(v)
    n = numel(categories(removecats(v)));
    return;
end
if isstring(v)
    v = v(~ismissing(v));
end
n = numel(unique(v));
end
