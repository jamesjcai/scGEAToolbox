function [t, srcname] = i_pickworkspacetable(parentfig, nrows, prompt)
% I_PICKWORKSPACETABLE - Pick a table variable from the base workspace.
%
%   [t, srcname] = gui.i_pickworkspacetable(parentfig)
%   [t, srcname] = gui.i_pickworkspacetable(parentfig, nrows, prompt)
%
% Inputs:
%   parentfig - parent figure for the dialogs (default: [])
%   nrows     - required number of rows, e.g. sce.NumCells; [] accepts any
%               number of rows (default: [])
%   prompt    - instruction shown as a wrapped label above the list; the
%               window title stays 'Workspace Table' either way (default:
%               explains the one-row-per-cell requirement)
%
% Outputs:
%   t       - the chosen table, [] when cancelled or when nothing qualified
%   srcname - name of the variable chosen, '' when cancelled
%
% Only variables of class table are offered, and when nrows is given only
% those with that many rows, so a table that cannot line up one row per cell
% never reaches the list. When nothing qualifies, saying so and why is more
% use than an empty list: the row counts that were found are named, since a
% table of the wrong length is the usual reason.
%
% see also: gui.i_readtablefile, gui.sc_cellattribeditor,
%           gui.callback_AssignCellTypeFromAttrib

if nargin < 1, parentfig = []; end
if nargin < 2, nrows = []; end
if nargin < 3 || isempty(prompt)
    if isempty(nrows)
        prompt = 'Select a table variable in the MATLAB workspace.';
    else
        prompt = sprintf(['Select a table variable in the MATLAB workspace. ', ...
            'It must have one row per cell (%d rows), in the same order as ', ...
            'the cells in this dataset.'], nrows);
    end
end

t = [];
srcname = '';

if ismcc || isdeployed
    gui.myWarndlg(parentfig, ['There is no base workspace to read from in ', ...
        'a deployed application. Load the table from a file instead.']);
    return;
end

a = evalin('base', 'whos');
if isempty(a)
    a = struct('name', {}, 'class', {}, 'size', {});
end
istbl = strcmp({a.class}, 'table');
names = {a(istbl).name};
if isempty(names)
    gui.myWarndlg(parentfig, ['There are no table variables in the base ', ...
        'workspace. Create one there, or load the table from a file ', ...
        'instead.']);
    return;
end

% Row counts come from whos, so a large table is not copied into this
% function merely to be measured.
nrowof = cellfun(@(s) s(1), {a(istbl).size});
if ~isempty(nrows)
    keep = nrowof == nrows;
    if ~any(keep)
        gui.myWarndlg(parentfig, sprintf(['None of the %d table variable(s) ', ...
            'in the base workspace has one row per cell (%d rows needed). ', ...
            'Row counts found: %s.'], numel(names), nrows, ...
            strjoin(string(unique(nrowof)), ', ')));
        return;
    end
    names = names(keep);
    nrowof = nrowof(keep);
end

% Naming the columns as well as the row count is what makes the list usable
% when the workspace holds several same-sized tables, which is the common case
% after a couple of exports.
items = cell(1, numel(names));
for k = 1:numel(names)
    cols = evalin('base', sprintf('%s.Properties.VariableNames', names{k}));
    collist = strjoin(string(cols), ', ');
    if strlength(collist) > 60
        collist = extractBefore(collist, 58) + '...';
    end
    items{k} = sprintf('%s  -  %d rows: %s', names{k}, nrowof(k), collist);
end

if gui.i_isuifig(parentfig)
    [indx, tf] = gui.myListdlg(parentfig, items, 'Workspace Table', [], ...
        false, true, [460, 220], prompt);
else
    [indx, tf] = listdlg('PromptString', {prompt}, 'SelectionMode', ...
        'single', 'ListString', items, 'ListSize', [400, 300]);
end
if tf ~= 1 || isempty(indx), return; end

srcname = names{indx};
try
    t = evalin('base', srcname);
catch ME
    gui.myErrordlg(parentfig, ME.message, ME.identifier);
    t = [];
    srcname = '';
    return;
end

% The variable could have been replaced between whos and now.
if ~istable(t) || (~isempty(nrows) && height(t) ~= nrows)
    gui.myErrordlg(parentfig, sprintf(['Variable "%s" is no longer a ', ...
        'table with the expected number of rows.'], srcname));
    t = [];
    srcname = '';
end
end
