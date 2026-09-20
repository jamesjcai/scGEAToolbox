function [listitems, wsvars, wsinfo] = i_classlistitems(sce, allowsingle)
%I_CLASSLISTITEMS  Categorical grouping variables offered by the class pickers.
%
%   [listitems, wsvars, wsinfo] = i_classlistitems(sce, allowsingle)
%
% Shared by GUI.I_SELECT1CLASS (pick one) and GUI.I_SELECTNCLASS (pick
% several and cross them). Kept in one place so the two dialogs can never
% offer different variables - they are meant to be the same list, and the
% only difference between the callers is how many items may be chosen.
%
% ALLOWSINGLE true (default) lists a variable whenever it is populated.
% False lists it only when it has more than one level, i.e. only when it
% could actually separate cells into groups.
%
% WSVARS / WSINFO are the base-workspace variables whose length matches
% sce.NumCells, as returned by WHOS and STRUCT2CELL respectively. They are
% non-empty only when 'Workspace Variable...' made it into LISTITEMS, and
% the caller needs them to resolve that choice.

if nargin < 2 || isempty(allowsingle), allowsingle = true; end

wsvars = [];
wsinfo = {};

listitems = {'Current Class (C)'};
listitems = i_add(listitems, sce.c_cluster_id,    'Cluster ID',       allowsingle);
listitems = i_add(listitems, sce.c_cell_type_tx,  'Cell Type',        allowsingle);
listitems = i_add(listitems, sce.c_cell_cycle_tx, 'Cell Cycle Phase', allowsingle);
listitems = i_add(listitems, sce.c_batch_id,      'Batch ID',         allowsingle);

a = evalin('base', 'whos');
b = struct2cell(a);
v = false(length(a), 1);
for k = 1:length(a)
    if max(a(k).size) == sce.NumCells && min(a(k).size) == 1
        v(k) = true;
    end
end
if any(v)
    wsvars = a(v);
    wsinfo = b(:, v);
    listitems = [listitems, 'Workspace Variable...'];
end

% Named cell attributes, when they are labels rather than measurements, go
% last and under a divider. An imported Seurat or RDS object puts every
% meta.data column it does not recognize into LIST_CELL_ATTRIBUTES, so this
% is where a categorical annotation ends up when it is not one of the
% fields above, and the two groups are worth telling apart on sight.
%
% The resolvers read these back by name through SCE.GETCELLATTRIBUTE, so an
% attribute sharing a name with an item already listed is skipped: it could
% never be reached past the switch that matches the built-in first.
attribitems = i_attributeitems(sce, listitems, allowsingle);
if ~isempty(attribitems)
    listitems = [listitems, {i_classlistdivider()}, attribitems];
end

end


function items = i_attributeitems(sce, taken, allowsingle)
% The names of every cell attribute that qualifies as a grouping variable,
% in the order they are stored. LIST_CELL_ATTRIBUTES is a flat
% {name1, value1, name2, value2, ...} cell.
items = {};
names = sce.list_cell_attributes(1:2:end);
for k = 1:numel(names)
    name = char(string(names{k}));
    if isempty(name) || any(strcmp(name, taken)) || any(strcmp(name, items))
        continue;
    end
    if 2*k > numel(sce.list_cell_attributes), continue; end   % name with no value
    if ~pkg.i_isgroupingvar(sce.list_cell_attributes{2*k}, sce.NumCells, allowsingle)
        continue;
    end
    items = [items, {name}]; %#ok<AGROW>
end
end


function listitems = i_add(listitems, attrib, name, allowsingle)
% One entry, added when the attribute is populated - and, when single-level
% variables are not wanted, only when it has something to separate.
if isempty(attrib), return; end
if ~allowsingle && numel(unique(attrib)) <= 1, return; end
listitems = [listitems, name];
end
