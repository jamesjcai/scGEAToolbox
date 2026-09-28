function [ctypelist, matched, primarytypes, resolved] = i_subtypecandidates(sce)
%I_SUBTYPECANDIDATES Primary cell types in a dataset that have subtype markers.
%
%   [ctypelist, matched, primarytypes, resolved] = pkg.i_subtypecandidates(sce)
%
% Outputs:
%   ctypelist    - the primary types present in sce.c_cell_type_tx; empty
%                  when the data has none of them
%   matched      - per cell, the primary type its label belongs to, "" for none
%   primarytypes - every primary type the marker table knows about
%   resolved     - per cell, true when the label already names a subtype of
%                  that primary, so a run will re-derive it rather than
%                  subdivide a cell that had no subtype yet
%
% One answer for both the menu and the dialog: GUI.I_UPDATEANNOTATEMENU asks
% whether there is anything to annotate before enabling the menu item, and
% GUI.CALLBACK_SUBTYPEANNOTATION asks the same question again to build its list.
% The marker table's cell type column is cached because the menu asks every
% time the Annotate menu is opened.
%
% RESOLVED is for reporting, not for filtering. A dataset annotated only as
% "Plasma cells" reaches the primary type "B cells" through
% PKG.I_SUBTYPEOVERLAP, and that run is worth offering: those cells are B
% cells whose subtype was decided by a competition among primary types rather
% than among B cell subtypes, and re-deriving it is the point. What RESOLVED
% buys is a dialog that can say so before the user commits to it.
%
% see also: pkg.i_matchprimarytype, pkg.i_subtypeoverlap,
%           gui.callback_SubtypeAnnotation, gui.i_updateannotatemenu

ctypelist = strings(0, 1);
matched = strings(0, 1);
resolved = false(0, 1);
primarytypes = in_primarytypes();

if nargin < 1 || isempty(sce) || ~isa(sce, 'SingleCellExperiment'), return; end
if sce.NumCells == 0 || isempty(sce.c_cell_type_tx), return; end
if isempty(primarytypes), return; end

matched = pkg.i_matchprimarytype(sce.c_cell_type_tx, primarytypes);

% Already subtyped counts only against the primary the label itself names:
% "Plasma cells" is a resolved B cell, and says nothing about the T cells in
% the same dataset.
[oprimary, ~, isoverlap] = pkg.i_subtypeoverlap(sce.c_cell_type_tx);
resolved = isoverlap & oprimary == matched;

ctypelist = unique(matched(strlength(matched) > 0));
end

function primarytypes = in_primarytypes()
% The CellType column of the bundled subtype marker table. Cached: reading the
% spreadsheet on every menu open would be felt. The cache is dropped when the
% file changes on disk, so adding a cell type to it takes effect without
% restarting MATLAB.

persistent cached cachedstamp

primarytypes = strings(0, 1);
pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '..', 'assets', 'PanglaoDB', 'cellsubtypes.xlsx');
if ~exist(pth, 'file')
    pth = which('cellsubtypes.xlsx');   % deployed layout
end
if isempty(pth) || ~exist(pth, 'file'), return; end

d = dir(pth);
stamp = [d.datenum, d.bytes];
if ~isempty(cached) && isequal(stamp, cachedstamp)
    primarytypes = cached;
    return;
end

try
    t = readtable(pth, 'TextType', 'string');
catch
    % An unreadable marker table leaves the menu item disabled rather than
    % failing on menu open.
    return;
end
if ~ismember('CellType', t.Properties.VariableNames), return; end

primarytypes = unique(t.CellType(strlength(t.CellType) > 0));
cached = primarytypes;
cachedstamp = stamp;
end
