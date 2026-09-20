function sortby = i_askgrouporder(parentfig, defaultby)
%I_ASKGROUPORDER Ask which order to list or draw a grouping variable in.
%
%   sortby = gui.i_askgrouporder(parentfig)
%   sortby = gui.i_askgrouporder(parentfig, defaultby)
%
%   Returns "count", "natural" or "none", ready for
%   GUI.I_ORDERGROUPLEVELS, or "" when the user cancelled.
%
%   defaultby is the button offered first, given the same way. It defaults
%   to "count": a long annotation is worked through by size, and the
%   groups worth looking at are the populous ones.
%
%   One dialog for every place the user picks a group order, so the wording
%   is the same whether the groups are being chosen for an analysis or
%   re-sorted on a plot that is already up.
%
%   See also GUI.I_ORDERGROUPLEVELS, GUI.I_SELECTGROUPSUBSET.

if nargin < 2 || isempty(defaultby), defaultby = "count"; end
if nargin < 1, parentfig = []; end

labels = {'Number of cells', 'Alphabetic', 'Unsorted'};
switch string(defaultby)
    case "natural"
        defaultlabel = labels{2};
    case "none"
        defaultlabel = labels{3};
    otherwise
        defaultlabel = labels{1};
end

switch gui.myQuestdlg(parentfig, 'Order groups by:', '', labels, defaultlabel)
    case 'Number of cells'
        sortby = "count";
    case 'Alphabetic'
        sortby = "natural";
    case 'Unsorted'
        sortby = "none";
    otherwise
        sortby = "";
end

end
