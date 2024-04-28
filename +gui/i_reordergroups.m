function [c, cL, noanswer, newidx] = i_reordergroups(thisc, preorderedcL, parentfig)

if nargin < 3, parentfig = []; end
if nargin < 2, preorderedcL = []; end
noanswer = true;
[c, cL] = grp2idx(thisc);
newidx=1:numel(cL);
if numel(cL) == 1
    %errordlg('Only one cell type or cluster.');
    noanswer = false;
    return;
end

[answer] = questdlg('Manually order groups?', '', ...
    'Yes', 'No', 'Cancel', 'No');
if isempty(answer), return; end
    switch answer
        case 'Yes'
            if ~isempty(preorderedcL)
                [newidx] = gui.i_selmultidlg(cL, preorderedcL, parentfig);
            else
                [newidx] = gui.i_selmultidlg(cL, natsort(cL), parentfig);
            end
            if length(newidx) ~= length(cL)
                noanswer = true;
                warndlg('Please select all items.','','modal');
                return;
            end
            cx = c;
            for k = 1:length(newidx)
                c(cx == newidx(k)) = k;
            end
            cL = cL(newidx);
            noanswer = false;
        case 'No'
            noanswer = false;
    end
end
