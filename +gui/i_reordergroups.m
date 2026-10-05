function [c, cL, noanswer, newidx] = i_reordergroups(thisc, preorderedcL, ...
    parentfig)

if nargin < 3, parentfig = []; end
if nargin < 2, preorderedcL = []; end
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end
noanswer = true;
[c, cL] = findgroups(string(thisc));
newidx = 1:numel(cL);
if isscalar(cL)
    noanswer = false;
    return;
end

[answer] = gui.myQuestdlg(parentfig, 'Manually order groups?', '', ...
{'Yes', 'No', 'Cancel'}, 'No');
if isempty(answer), return; end
switch answer
    case 'Yes'
        if isempty(preorderedcL)
            preorderedcL = natsort(cL);
        end
        [newidx] = gui.i_selmultidialog(cL, preorderedcL, parentfig);
        if length(newidx) ~= length(cL)
            gui.myWarndlg(parentfig, 'Please select all items.');
            return;
        end
    case 'No'
        % No manual order is not FINDGROUPS order: that is plain character
        % order, "Cluster 10" before "Cluster 2", where every group list
        % the user has been shown is natural-sorted. PREORDEREDCL, when
        % given, is the order a list was shown in; levels it leaves out
        % follow in FINDGROUPS order.
        if isempty(preorderedcL)
            [~, newidx] = natsort(cL);
        else
            [~, newidx] = ismember(string(preorderedcL), cL);
            newidx = newidx(newidx > 0);
            newidx = [newidx(:); setdiff((1:numel(cL)).', newidx(:))];
        end
    otherwise
        % 'Cancel'
        return;
end

newidx = newidx(:).';
cx = c;
for k = 1:length(newidx)
    c(cx == newidx(k)) = k;
end
cL = cL(newidx);
noanswer = false;
end
