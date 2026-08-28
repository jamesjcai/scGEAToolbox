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

end


function listitems = i_add(listitems, attrib, name, allowsingle)
% One entry, added when the attribute is populated - and, when single-level
% variables are not wanted, only when it has something to separate.
if isempty(attrib), return; end
if ~allowsingle && numel(unique(attrib)) <= 1, return; end
listitems = [listitems, name];
end
