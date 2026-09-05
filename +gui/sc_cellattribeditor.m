function [sce, needupdate] = sc_cellattribeditor(sce, addnew, parentfig)
if nargin<3, parentfig = []; end
if nargin<2, addnew = false; end
if ~isempty(parentfig)
    figure(parentfig);
    cleanupObj = onCleanup(@() figure(parentfig));
end

needupdate = false;


if ~addnew    % edit
    baselistitems = {'Cell Cycle Phase', ...
        'Cell Type', 'Cluster ID', ...
        'Batch ID', 'Cell ID'};
    listitems = [baselistitems, ...
        sce.list_cell_attributes(1:2:end)];
    listitems = listitems(~cellfun(@isempty, listitems));

        if gui.i_isuifig(parentfig)
            [indx2, tf2] = gui.myListdlg(parentfig, listitems, 'Select Cell Attribute:');
        else
            [indx2, tf2] = listdlg('PromptString', ...
                {'Select Cell Attribute:'}, ...
                'SelectionMode', 'single', 'ListString', listitems, ...
                'ListSize', [220, 300]);
        end

    if tf2 == 1
        clabel = listitems{indx2};
        thisc = [];
        switch clabel
            case 'Cluster ID'   % cluster id
                thisc = sce.c_cluster_id;
            case 'Batch ID'       % batch id
                thisc = sce.c_batch_id;
            case 'Cell Type'  % cell type
                thisc = sce.c_cell_type_tx;
            case 'Cell Cycle Phase' % cell cycle
                thisc = sce.c_cell_cycle_tx;
            case 'Cell ID'
                thisc = sce.c_cell_id;
            otherwise
                [y,idx] = ismember(clabel, sce.list_cell_attributes(1:2:end));
                if y
                    thisc = sce.list_cell_attributes{2*idx};
                end
        end
    else
        return;
    end

    % Values can be edited one distinct value at a time, typed in one line per
    % cell, or taken from a column of a per-cell table - this dataset's own cell
    % attribute table, a delimited text file, or a table variable in the base
    % workspace, the same sources gui.callback_AssignCellTypeFromAttrib offers.
    % Distinct-value editing is the only one that does not grow with the number
    % of cells, so it leads for a long attribute.
    canrelabel = in_canrelabel(thisc);
    islong = numel(thisc) > in_maxrawvalues();
    tsce = in_scetable(sce, clabel);
    srctag = in_picksource(parentfig, sprintf( ...
        'Where should the values for "%s" come from?', clabel), ...
        canrelabel, islong, 'Cell Attribute Values', width(tsce) > 0);
    switch srctag
        case 'relabel'
            [sce, needupdate] = in_relabelattrib(sce, clabel, thisc, parentfig);
            return;
        case {'sce', 'file', 'workspace'}
            if ~in_confirmnewannotation(sce, clabel, parentfig), return; end
            [sce, needupdate] = in_replaceattribfromtable(sce, clabel, thisc, ...
                parentfig, srctag, tsce);
            return;
        case 'manual'
            if ~in_confirmnewannotation(sce, clabel, parentfig), return; end
            if ~in_confirmrawedit(parentfig, numel(thisc)), return; end
        otherwise
            return;
    end

%    if ~strcmp('Yes', gui.myQuestdlg(parentfig, ...
%            'It may take a while to load values. Continue?'))
%        return;
%    end
    % tic;
    if gui.i_isuifig(parentfig)
        %        x = gui.myInputdlg({sprintf('Attribute Name: %s\n%s',clabel, 'Attribute Values:')}, ...
        %                          'Attribute Editor', {char(string(thisc))}, parentfig);

        % assignin("base","thisc",thisc);
        % assignin("base","clabel",clabel);

        x = gui.myTextareadlg(parentfig, {'Attribute Name', 'Attribute Values'},...
                      'Attribute Editor', {clabel, string(thisc)}, [false, true]);
        % assignin("base","x",x);
        if ~isempty(x)
            x(1)=[];
        end
    else
        x = inputdlg(sprintf('Attribute Name: %s\n%s',clabel, 'Attribute Values:'), ...
                          'Attribute Editor', [15 80], {char(string(thisc))});

    end
    % toc;

else    % add new

    % Values can be typed in or taken from a column of a per-cell table: this
    % dataset's own cell attribute table, a delimited text file, or a table
    % variable in the base workspace.
    tsce = in_scetable(sce);
    srctag = in_picksource(parentfig, ...
        'Where should the values for the new attribute come from?', ...
        false, false, 'New Cell Attribute', width(tsce) > 0);
    switch srctag
        case {'sce', 'file', 'workspace'}
            [sce, needupdate] = in_addattribsfromtable(sce, parentfig, srctag, tsce);
            return;
        case 'manual'
            if ~in_confirmrawedit(parentfig, sce.NumCells), return; end
        otherwise
            return;
    end

    if gui.i_isuifig(parentfig)
        x = gui.myTextareadlg(parentfig, {'Attribute Name', 'Attribute Values'},...
                      'Attribute Editor', {'New_Attribute', ("Value_"+string(1:sce.NumCells))'});
    else
        x = inputdlg({'Attribute Name', 'Attribute Values'},...
                      'Attribute Editor', [1 80; 15 80]);
    end

    %{
    if gui.i_isuifig(parentfig)
        % x = gui.myInputdlg({'Attribute Name','Attribute Values'},...
        %              'Attribute Editor', {''}, parentfig); % Assuming default is empty cell
        % x = gui.myTextareadlg(parentfig, '', 'Attribute Name'); % Assuming default is empty cell
        x = gui.myTextareadlg(parentfig, {'Attribute Name','Attribute Values'},...
                      'Attribute Editor', {'New_Attribute', ("Value_"+string(1:sce.NumCells))'});
    else
        x = inputdlg({'Attribute Name','Attribute Values'},...
                      'Attribute Editor', [1 80; 15 80]);
    end
    %}
    % {'new_attrib', char(string([1:sce.NumCells]'))});
end

if isempty(x), return; end

if addnew
    if isempty(x{1})
        gui.myWarndlg(parentfig, 'Attribute Name cannot be empty.');
        return;
    end
    if isempty(x{2})
        gui.myWarndlg(parentfig, 'Attribute Values cannot be empty.');
        return;
    end
else
    if isempty(x{1})       % when add new - x{1} is the values
        gui.myWarndlg(parentfig, 'Attribute Values cannot be empty.');
        return;
    end
end

answer = gui.myQuestdlg(parentfig, 'What is the data type of attribute values?', ...
	    'Data Type', ...
	    {'String','Numeric','Cancel'},'String');
switch answer
        case 'String'
            if addnew
                clabel = strtrim(x{1});
                newthisc = strtrim(string(trimBottomEmpty(x{2})));
            else
                newthisc = strtrim(string(trimBottomEmpty(x{1})));    % when add new - x{1} is the values
            end
        case 'Numeric'
            if addnew
                clabel = strtrim(x{1});
                newthisc = str2double(string(trimBottomEmpty(x{2})));
            else
                newthisc = str2double(string(trimBottomEmpty(x{1})));   % when add new - x{1} is the values
            end
        case 'Cancel'
            return;
end


if addnew
    if size(newthisc,1) ~= sce.NumCells
        gui.myWarndlg(parentfig, ...
            'Attribute length is not equal to the number of cells.');
        return;
    end
else
    if ~isequal(size(newthisc), size(thisc))
        gui.myWarndlg(parentfig, 'Attribute length changed.');
        return;
    end
end

if addnew
    clabel = matlab.lang.makeValidName(clabel);
    existinglabels = sce.list_cell_attributes(1:2:end);
    if ismember(clabel, existinglabels)
        gui.myWarndlg(parentfig, 'Cell Attribute Name Existing.');
        return;
    end
    sce.list_cell_attributes = [sce.list_cell_attributes, ...
        {clabel, newthisc(:)}];
    gui.myHelpdlg(parentfig, 'Cell Attribute Added.');
    needupdate = true;
else
    [sce, stashname] = in_setattribvalue(sce, clabel, newthisc);
    gui.myHelpdlg(parentfig, ...
        "Cell Attribute Changed." + gui.i_stashnotice(stashname));
    needupdate = true;
end
end


function src = in_picksource(parentfig, prompt, allowrelabel, preferrelabel, ...
        dlgtitle, allowsce)
% Ask where the attribute values come from. '' when cancelled.
%
% allowrelabel adds distinct-value editing to the choices, and preferrelabel
% makes it the pre-selected one, which is what a long attribute wants. allowsce
% adds this dataset's own cell attribute table, and is false when that table has
% no column to offer - every attribute is the wrong length, or the only one
% there is the attribute being edited.
%
% The instruction goes in GUI.MYLISTDLG's prompt label, which wraps, rather
% than in its Title, which the window title bar clips without saying so. Each
% choice still names what it reads and in what shape, so the list stands on
% its own if the prompt is skimmed. The dialog is sized to the text rather
% than left at the tall, narrow default meant for long lists.

if nargin < 5 || isempty(dlgtitle), dlgtitle = 'Cell Attribute Values'; end
if nargin < 6, allowsce = false; end

src = '';
relabeltag = 'Rename the distinct values (edit labels one by one)';
typetag = 'Type in the values (one line per cell)';
scetag = 'Select a column from the cell attribute table';
loadtag = 'Read a column from a table file (CSV, TSV or TXT)';
wstag = 'Read a column from a table variable in the workspace';

% The dataset's own attribute table goes directly after typing the values in:
% both take what is already in this session, ahead of the two choices that
% bring a table in from outside it.
items = {typetag, loadtag, wstag};
if allowsce
    items = {typetag, scetag, loadtag, wstag};
end
prefersel = typetag;
if allowrelabel
    items = [{relabeltag}, items];
    if preferrelabel, prefersel = relabeltag; end
end

if gui.i_isuifig(parentfig)
    [indx, tf] = gui.myListdlg(parentfig, items, dlgtitle, prefersel, false, ...
        true, [480, 200], prompt);
else
    [indx, tf] = listdlg('PromptString', {prompt}, 'SelectionMode', 'single', ...
        'ListString', items, 'ListSize', [460, 100], ...
        'InitialValue', find(strcmp(items, prefersel), 1));
end
if tf ~= 1 || isempty(indx), return; end

switch items{indx}
    case relabeltag
        src = 'relabel';
    case scetag
        src = 'sce';
    case loadtag
        src = 'file';
    case wstag
        src = 'workspace';
    otherwise
        src = 'manual';
end
end

function [t, srcname] = in_gettable(srctag, parentfig, nrows, tsce)
% Fetch a per-cell table from whichever source in_picksource returned, so the
% two importers below do not each need to know about all of them. tsce is this
% dataset's own attribute table, already built by the caller to decide whether
% to offer it at all.

switch srctag
    case 'sce'
        t = tsce;
        srcname = 'the cell attribute table';
    case 'workspace'
        [t, srcname] = gui.i_pickworkspacetable(parentfig, nrows);
    otherwise
        [t, srcname] = gui.i_readtablefile(parentfig, nrows);
end
end

function t = in_scetable(sce, excludelabel)
% Build a per-cell table out of this dataset's own attributes: the standard
% fields followed by everything in list_cell_attributes, the same content the
% View/Export Cell Attribute Table shows.
%
% Unlike PKG.I_MAKEATTRIBUTESTABLE, each column keeps the class it is stored in
% rather than being converted to text, so copying one numeric attribute onto
% another does not send its values through str2double on the way.
%
% excludelabel drops one attribute by its display name, which is how the
% attribute being edited is kept from being offered as a source for itself.

if nargin < 2, excludelabel = ''; end

t = table();
n = sce.NumCells;
if n == 0, return; end

names = [{'Cell ID', 'Batch ID', 'Cluster ID', 'Cell Type', ...
    'Cell Cycle Phase'}, sce.list_cell_attributes(1:2:end)];
vals = [{sce.c_cell_id, sce.c_batch_id, sce.c_cluster_id, ...
    sce.c_cell_type_tx, sce.c_cell_cycle_tx}, ...
    sce.list_cell_attributes(2:2:end)];

keep = false(1, numel(names));
for k = 1:numel(names)
    if strcmp(names{k}, excludelabel), continue; end
    v = vals{k};
    if ischar(v)
        v = string(cellstr(v));   % char matrix, as older .mat files store cell ids
    elseif iscell(v)
        v = string(v);
    end
    if numel(v) ~= n, continue; end
    vals{k} = v(:);
    keep(k) = true;
end
if ~any(keep), return; end

% Table variable names have to be valid identifiers, and the standard fields
% carry spaces. IN_REPLACEATTRIBFROMTABLE ignores spaces and underscores when it
% pre-selects a column, so 'Cell Type' still matches CellType here.
colnames = matlab.lang.makeUniqueStrings( ...
    matlab.lang.makeValidName(names(keep)));
vals = vals(keep);
t = table(vals{:}, 'VariableNames', colnames);
end

function n = in_maxdistinct()
% More distinct values than this and the relabel dialog is itself a long list to
% work through, so the attribute is not worth offering it for.

n = 1000;
end

function n = in_maxrawvalues()
% Above this many values the one-line-per-cell editor is too slow to open.
% uitextarea load time grows faster than the number of lines it is given: 20k
% values appear at once, 100k take about half a minute, 200k several minutes.

n = 20000;
end

function tf = in_canrelabel(thisc)
% Distinct-value editing only means something for a grouping variable. A column
% of per-cell IDs has one distinct value per cell and nothing to relabel, and a
% continuous measurement has too many rows to show.

tf = false;
if isempty(thisc), return; end

if isnumeric(thisc)
    uniqvals = unique(thisc(:));
else
    uniqvals = unique(string(thisc));
end

nuniq = numel(uniqvals);
if nuniq == 0 || nuniq > in_maxdistinct() || nuniq >= numel(thisc), return; end

if isnumeric(thisc)
    % The dialog groups and returns values as text, so a numeric attribute can
    % only go through it when its values survive the round trip. Cluster ids do;
    % a measurement that takes more digits than string() prints does not. This
    % is checked after the count, which rules out continuous data far cheaper.
    if ~isequaln(str2double(string(uniqvals)), uniqvals), return; end
end

tf = true;
end

function tf = in_confirmrawedit(parentfig, nval)
% Warn before opening the one-line-per-cell editor on more values than it can
% load promptly. It is still allowed: the wait is the only cost.

tf = true;
if nval <= in_maxrawvalues(), return; end

answer = gui.myQuestdlg(parentfig, sprintf(['Loading %d values into the editor ' ...
    'can take minutes. Loading a table from a file instead is immediate. ' ...
    'Continue anyway?'], nval), 'Attribute Editor', {'Yes', 'No'}, 'No');
tf = strcmp(answer, 'Yes');
end

function tf = in_confirmnewannotation(sce, clabel, parentfig)
% Warn before every cell type value is replaced at once, the same way the
% annotation callbacks do. Only the wholesale paths ask: relabeling distinct
% values renames the annotation that is already there rather than supplying a
% different one, and asking on every rename would be noise. The stash in
% IN_SETATTRIBVALUE happens either way, so nothing is lost on a rename.

tf = true;
if ~strcmp(clabel, 'Cell Type'), return; end
tf = gui.i_confirmoverwritecelltype(parentfig, sce);
end

function [sce, needupdate] = in_relabelattrib(sce, clabel, thisc, parentfig)
% Edit the attribute one distinct value at a time. Values that parse as numbers
% are written back as numbers when that is what the attribute held, so relabeling
% a numeric cluster id 1,2,3 to 3,2,1 does not turn it into text; renaming those
% same ids to names does, which is the point of the edit.
%
% Cells that were NaN come back from the dialog as <missing>, which str2double
% turns into NaN again. That is the value they already had, so it does not count
% as a failure to parse - without the ismissing test below, one NaN anywhere in
% a numeric attribute would force the whole attribute to text on any edit.

needupdate = false;
newvals = gui.i_relabeldlg(parentfig, thisc, ...
    sprintf('Edit Distinct Values of "%s"', clabel));
if isempty(newvals), return; end

if isnumeric(thisc)
    numvals = str2double(newvals);
    if ~any(isnan(numvals) & ~ismissing(newvals))
        newvals = numvals;
    end
end

[sce, stashname] = in_setattribvalue(sce, clabel, newvals);
needupdate = true;

gui.myHelpdlg(parentfig, ...
    "Cell Attribute Changed." + gui.i_stashnotice(stashname));
end

function [sce, needupdate] = in_replaceattribfromtable(sce, clabel, thisc, ...
        parentfig, srctag, tsce)
% Replace the values of an existing attribute with a column of a per-cell table,
% read either from a file or from a table variable in the base workspace.
%
% The column is coerced to the class of the values it replaces rather than taken
% as the table holds it: Cluster ID has to stay numeric and Cell Type has to stay
% text, whatever the table happens to carry. An attribute that is empty so far has
% no class to match, so there the table's own type is kept.

needupdate = false;
[t, srcname] = in_gettable(srctag, parentfig, sce.NumCells, tsce);
if isempty(t), return; end
if width(t) == 0
    gui.myWarndlg(parentfig, sprintf('There are no columns in %s to choose from.', srcname));
    return;
end

% Pre-select a column named after the attribute: the usual case is a table
% exported from this dataset, edited elsewhere, and brought back. Spaces and
% underscores are ignored in the match, since a header like "Cell Type" reaches
% here as CellType once readtable has made it a valid identifier.
colnames = t.Properties.VariableNames;
squash = @(v) lower(regexprep(string(v), '[^a-zA-Z0-9]', ''));
prefersel = find(squash(colnames) == squash(clabel), 1);
[indx, tf] = gui.myListdlg(parentfig, colnames, ...
    sprintf('Select Column of Values for "%s":', clabel), prefersel, false);
if tf ~= 1 || isempty(indx), return; end

v = t.(colnames{indx});
if isempty(thisc)
    if ~isnumeric(v), v = string(v); end
elseif isnumeric(thisc)
    if ~isnumeric(v)
        v = str2double(string(v));
        if all(isnan(v))
            gui.myErrordlg(parentfig, sprintf(['Column "%s" holds text, but %s ' ...
                'takes numeric values.'], colnames{indx}, clabel));
            return;
        end
    end
else
    v = string(v);
end

% Match the orientation of what is being replaced; the table column is a column
% vector but not every per-cell field is stored as one.
v = v(:);
if isrow(thisc) && numel(thisc) > 1
    v = v.';
end

[sce, stashname] = in_setattribvalue(sce, clabel, v);
needupdate = true;

msg = sprintf('%s replaced with column "%s" of %s.', ...
    clabel, colnames{indx}, srcname);
gui.myHelpdlg(parentfig, msg + gui.i_stashnotice(stashname));
end

function [sce, stashname] = in_setattribvalue(sce, clabel, newthisc)
% Write values to the attribute named by clabel, whether it is one of the
% standard per-cell fields or an entry of list_cell_attributes. Returns the
% name of the 'old_cell_type_N' attribute the previous cell type labels were
% stashed in, or "" when nothing was stashed, so the caller can say where
% they went (see GUI.I_STASHNOTICE).

stashname = "";
switch clabel
    case 'Cluster ID'
        sce.c_cluster_id = newthisc;
    case 'Batch ID'
        sce.c_batch_id = newthisc;
    case 'Cell Type'
        % Every path through the editor lands here, so one stash covers the
        % relabel, typed-in and from-file edits alike (see
        % PKG.I_STASHCELLTYPEHISTORY).
        stashname = pkg.i_stashcelltypehistory(sce);
        sce.c_cell_type_tx = newthisc;
    case 'Cell Cycle Phase'
        sce.c_cell_cycle_tx = newthisc;
    case 'Cell ID'
        sce.c_cell_id = newthisc;
    otherwise
        if sce.hasCellAttribute(clabel)
            sce.setCellAttribute(clabel, newthisc);
        end
end
end

function [sce, needupdate] = in_addattribsfromtable(sce, parentfig, srctag, tsce)
% Add one or more columns of a per-cell table as new cell attributes, read
% either from a file or from a table variable in the base workspace.
%
% The table carries its own column types, so unlike the typed-in path there is no
% data type question to answer: numeric columns are stored as numbers and
% everything else as text. More than one column can be added at a time, which is
% the point of importing a table rather than pasting one column of values.

needupdate = false;
[t, srcname] = in_gettable(srctag, parentfig, sce.NumCells, tsce);
if isempty(t), return; end
if width(t) == 0
    gui.myWarndlg(parentfig, sprintf('There are no columns in %s to add.', srcname));
    return;
end

colnames = t.Properties.VariableNames;
[indx, tf] = gui.myListdlg(parentfig, colnames, ...
    'Select Column(s) to Add as Cell Attributes:', [], true);
if tf ~= 1 || isempty(indx), return; end

% A column whose name is already taken is added under a suffixed name instead of
% overwriting the attribute that is already there.
existinglabels = sce.list_cell_attributes(1:2:end);
newlabels = matlab.lang.makeUniqueStrings( ...
    matlab.lang.makeValidName(colnames(indx)), existinglabels);

for k = 1:numel(indx)
    v = t.(colnames{indx(k)});
    if ~isnumeric(v)
        v = string(v);
    end
    sce.setCellAttribute(newlabels{k}, v(:));
end
needupdate = true;

gui.myHelpdlg(parentfig, sprintf('%d cell attribute(s) added from %s: %s.', ...
    numel(indx), srcname, strjoin(string(newlabels), ', ')));
end

function trimmed_text = trimBottomEmpty(input_text)
% TRIMBOTTOMEMPTY Trims whitespace from each line and removes empty lines from bottom
%
% SYNTAX:
%   trimmed_text = trimBottomEmpty(input_text)
%
% INPUT:
%   input_text - char array or string, can be multiline
%
% OUTPUT:
%   trimmed_text - char array with whitespace trimmed from each line
%                  and empty lines removed from bottom only
%
% EXAMPLE:
%   text = ['111'; '222'; '333'; '   '];
%   result = trimBottomEmpty(text);

lines = cellstr(input_text);     % Convert to cell array
lines = strtrim(lines);          % Trim each line

% Find the last non-empty line
last_non_empty = find(~cellfun('isempty', lines), 1, 'last');

if ~isempty(last_non_empty)
    lines = lines(1:last_non_empty);  % Keep only up to last non-empty line
else
    lines = {};  % All lines were empty
end

trimmed_text = char(lines);      % Convert back to char array
end
