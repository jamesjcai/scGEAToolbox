function callback_Violinplot(src, ~)
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
[thisc, ~] = gui.i_select1class(sce);
if isempty(thisc), return; end

[c, cL, noanswer] = gui.i_reordergroups(thisc);
if noanswer, return; end

% [c, cL] = grp2idx(thisc);
% [answer] = questdlg('Manually order groups?', '', ...
%     'Yes', 'No', 'Cancel', 'No');
% if isempty(answer), return; end
% switch answer
%     case 'Yes'
%         [newidx] = gui.i_selmultidlg(cL, natsort(cL));
%         if length(newidx) ~= length(cL)
%             warndlg('Please select all group items.', '');
%             return;
%         end
%         cx = c;
%         for k = 1:length(newidx)
%             c(cx == newidx(k)) = k;
%         end
%         cL = cL(newidx);
%     case 'No'
% 
%     case 'Cancel'
%         return;
%     otherwise
%         return;
% end

[glist] = gui.i_selectngenes(sce);
if isempty(glist)
    helpdlg('No gene selected.', '');
    return;
end
[Xt] = gui.i_transformx(sce.X);
% glist=glist(end:-1:1);

try

%f=gui.i_violinmatrix(Xt,sce.g,c,cL,"PRLR");

for k = 1:length(glist)
    y = Xt(sce.g == glist(k), :);
    [f] = gui.i_violinplot(y, cL(c), glist(k), true, cL);
    % f=gui.i_violinplot(sce.X,sce.g,c,cL,glist);
    p = f.Position;
    %p(2)=0;
    p(4) = p(4) * 1.5;
    f.Position = p;
    f.Visible = "on";
end
    catch ME
        if exist('f','var') && ishandle(f)
            close(f);
        end
        errordlg(ME.message);
    end
end