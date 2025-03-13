function [i1, i2, cL1, cL2] = i_select2smplgrps(sce, needreorder, parentfig)


if nargin < 3, parentfig = []; end
if nargin < 2, needreorder = true; end
i1 = 0;
i2 = 0;
cL1 = [];
cL2 = [];



if isa(sce, 'SingleCellExperiment')
    [thisc, ~] = gui.i_select1class(sce,[],[],[],parentfig);
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

listitems = string(cLi);
n = length(listitems);
if n < 2
    gui.myErrordlg(parentfig, 'Need at least two samples.');
    return;
end




% [ci, cLi] = grp2idx(thisc);
 if n == 2
     answer = 'Two Samples';
 elseif n > 2
     [answer] = gui.myQuestdlg(parentfig, 'Select two samples or two sample groups?', '', ...
             {'Two Samples', 'Two Sample Groups', 'Cancel'}, 'Two Samples');
 end
 if isempty(answer), return; end
 switch answer
     case 'Two Samples'

        if gui.i_isuifig(parentfig)
            [indxx, tfx] = gui.myListdlg(parentfig, listitems, 'Select two samples:');
        else
            [indxx, tfx] = listdlg('PromptString', {'Select two samples:'}, ...
                'SelectionMode', 'multiple', ...
                'ListString', listitems, ...
                'InitialValue', [n - 1, n], 'ListSize', [220, 300]);
        end


        if tfx == 1
            if numel(indxx) ~= 2
                gui.myErrordlg(parentfig, 'Please select 2 samples only');
                return;
            end
            [y1, idx1] = ismember(listitems(indxx(1)), cLi);
            [y2, idx2] = ismember(listitems(indxx(2)), cLi);
            assert(y1 & y2);
            i1 = ci == idx1;
            i2 = ci == idx2;
            cL1 = cLi(idx1);
            cL2 = cLi(idx2);
        else
            return;
        end
     case 'Two Sample Groups'
         answer = gui.myQuestdlg(parentfig, 'Select samples in group 1?','');
         if ~strcmp(answer, 'Yes'), return; end
         [newidx1] = gui.i_selmultidlg(cLi, [], parentfig);
         % if length(newidx1) == length(cLi)
         %     gui.myWarndlg(parentfig, 'Please select all group items.', '');
         %     return;
         % end
         cx = ci;
         ci = zeros(size(ci));
         for k = 1:length(newidx1)
             ci(cx == newidx1(k)) = 1;
         end
         newgrp1name = inputdlg('Input a name for group 1', '', ...
            [1, 50], {'Group1'});
            if ~isempty(newgrp1name)
                cL1 = newgrp1name;
            else
                return;
            end

         answer = gui.myQuestdlg(parentfig, 'Select samples in group 2?','');
         if ~strcmp(answer, 'Yes'), return; end
         [newidx2] = gui.i_selmultidlg(cLi, [], parentfig);
         for k = 1:length(newidx2)
             ci(cx == newidx2(k)) = 2;
         end
         newgrp2name = inputdlg('Input a name for group 2', '', ...
            [1, 50], {'Group2'});
            if ~isempty(newgrp2name)
                cL2 = newgrp2name;
            else
                return;
            end
         
         i1 = ci == 1;
         i2 = ci == 2;
         % cL1 = {'Group1'};
         % cL2 = {'Group2'};
     otherwise
         return;
 end
end