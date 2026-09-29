function lines = i_menushortcuts(fig)
%I_MENUSHORTCUTS One line per menu item with a keyboard shortcut.
%
%   lines = gui.i_menushortcuts(fig) returns a string column such as
%   "Ctrl+K   Cluster > Cluster Cells", sorted by key, for every uimenu
%   under FIG whose Accelerator is set. Empty when FIG has none.
%
%   See also GUI.CALLBACK_SHORTCUTSGUIDE.

lines = strings(0, 1);
if isempty(fig) || ~isgraphics(fig), return; end
items = findall(fig, 'Type', 'uimenu');
items = items(~cellfun(@isempty, {items.Accelerator}));
if isempty(items), return; end

keys = strings(numel(items), 1);
paths = strings(numel(items), 1);
for k = 1:numel(items)
    keys(k) = upper(string(items(k).Accelerator));
    paths(k) = in_menupath(items(k));
end
[keys, order] = sort(keys);
lines = "Ctrl+" + keys + "   " + paths(order);
end

function p = in_menupath(m)
% "Cluster > Cluster Cells": the item's label and its parents', with the
% mnemonic ampersands and the trailing ellipsis taken off.
parts = strings(0, 1);
while isa(m, 'matlab.ui.container.Menu')
    parts = [in_clean(m.Text); parts]; %#ok<AGROW>
    m = m.Parent;
end
p = strjoin(parts, " > ");
end

function t = in_clean(t)
t = string(t);
t = replace(t, "&&", char(0));
t = erase(t, "&");
t = replace(t, char(0), "&");
if endsWith(t, "...")
    t = extractBefore(t, strlength(t) - 2);
end
t = strip(t);
end
