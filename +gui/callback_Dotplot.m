function callback_Dotplot(src, ~)
    if isa(src, 'matlab.apps.AppBase')    
        [FigureHandle, sce] = xui.gui_getfigsce(src);
    else
        [FigureHandle, sce] = gui.gui_getfigsce(src);
    end
    [thisc, ~] = gui.i_select1class(sce,[],[],[],FigureHandle);
    if isempty(thisc), return; end
    
    [c, cL, noanswer] = gui.i_reordergroups(thisc, [], FigureHandle);
    if noanswer, return; end
    % [c, cL] = grp2idx(thisc);
    % 
    % [answer] = gui.myQuestdlg(FigureHandle, 'Manually order groups?', '', ...
    %     'Yes', 'No', 'Cancel', 'No');
    % if isempty(answer), return; end
    % switch answer
    %     case 'Yes'
    %         [newidx] = gui.i_selmultidlg(cL, natsort(cL));
    %         if length(newidx) ~= length(cL)
    %             gui.myWarndlg(FigureHandle, 'Please select all group items.', '');
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
    
    [glist] = gui.i_selectngenes(sce, [], FigureHandle);
    if isempty(glist)
        gui.myHelpdlg(FigureHandle, 'No gene selected.', '');
        return;
    end
    [Xt] = gui.i_transformx(sce.X, true, 3, FigureHandle);
    if isempty(Xt), return; end
    glist = glist(end:-1:1);
    
    % if length(glist) > 50
    % 
    %     answer = gui.myQuestdlg(FigureHandle, 'Output to PowerPoint?');
    %     switch answer
    %         case 'Yes'
    %             needpptx = true;
    %         case 'No'
    %             needpptx = false;
    %         otherwise
    %             return;
    %     end
    %     images = {};
    %     for k = 1:50:length(glist)
    %         k2 = min([length(glist), k + 50 - 1]);
    %         hFig = gui.i_dotplot(Xt, sce.g, c, cL, glist(k:k2), true);
    %         screensize = get(groot, 'Screensize');
    %         p = hFig.Position;
    %         p(2) = 0;
    %         p(4) = screensize(4) - 150;
    %         hFig.Position = p;
    %         if needpptx
    %             img1 = [tempname, '.png'];
    %             images = [images, {img1}];
    %             saveas(hFig, img1);
    %         end
    %     end
    %     if needpptx, gui.i_save2pptx(images); end
    % 
    % else
        try
            gui.i_dotplot(Xt, sce.g, c, cL, glist, true, 'Dotplot', FigureHandle);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        end
    % end
end