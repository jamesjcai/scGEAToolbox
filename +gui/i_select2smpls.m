function [i1, i2, cL1, cL2] = i_select2smpls(sce, needreorder, parentfig)
if nargin < 3, parentfig = []; end
if nargin < 2, needreorder = true; end
i1 = 0;
i2 = 0;
cL1 = [];
cL2 = [];



if isa(sce, 'SingleCellExperiment')
    [thisc, ~] = gui.i_select1class(sce);
    if isempty(thisc)
        %gui.myErrordlg(parentfig, 'Undefined');
        return;
    end
else % assume that sce input is thisc
    thisc = sce;
end

if isscalar(unique(thisc))
    gui.myWarndlg(parentfig, "Cannot compare with an unique group");
    return;
end

if needreorder
    [ci, cLi, noanswer] = gui.i_reordergroups(thisc, [], parentfig);
    if noanswer, return; end
else
    [ci, cLi] = grp2idx(thisc);
end


% [ci, cLi] = grp2idx(thisc);
% 
% [answer] = gui.myQuestdlg(parentfig, 'Manually order groups?', '', ...
%     'Yes', 'No', 'Cancel', 'No');
% if isempty(answer), return; end
% switch answer
%     case 'Yes'
%         [newidx] = gui.i_selmultidlg(cLi, natsort(cLi));
%         if length(newidx) ~= length(cLi)
%             gui.myWarndlg(parentfig, 'Please select all group items.', '');
%             return;
%         end
%         cx = ci;
%         for k = 1:length(newidx)
%             ci(cx == newidx(k)) = k;
%         end
%         cLi = cLi(newidx);
%     case 'No'
% 
%     case 'Cancel'
%         return;
%     otherwise
%         return;
% end



listitems = string(cLi);
n = length(listitems);
if n < 2
    gui.myErrordlg(parentfig, 'Need at least two groups.');
    return;
end

        if gui.i_isuifig(parentfig)
            [indxx, tfx] = gui.myListdlg(parentfig, listitems, 'Select two groups:');
        else
            [indxx, tfx] = listdlg('PromptString', {'Select two groups:'}, ...
                'SelectionMode', 'multiple', ...
                'ListString', listitems, ...
                'InitialValue', [n - 1, n], 'ListSize', [220, 300]);
        end


if tfx == 1
    if numel(indxx) ~= 2
        gui.myErrordlg(parentfig, 'Please select 2 groups');
        return;
    end
    [y1, idx1] = ismember(listitems(indxx(1)), cLi);
    [y2, idx2] = ismember(listitems(indxx(2)), cLi);
    assert(y1 & y2);
    i1 = ci == idx1;
    i2 = ci == idx2;
    cL1 = cLi(idx1);
    cL2 = cLi(idx2);

    %i1=ismember(ci,indxx(1));
    %i2=ismember(ci,indxx(2));
    %cL1=cLi(indxx(1));
    %cL2=cLi(indxx(2));
else
    return;
end
end