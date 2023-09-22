function [requirerefresh, highlightindex, newclassidenty] = xcallback_AnnotateSubTypes(src, ~)


mfolder = fileparts(mfilename('fullpath'));
requirerefresh = false;
highlightindex = [];
newclassidenty = [];

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
thisc = gui.i_select1class(sce);
if isempty(thisc), return; end
[ci, cLi] = grp2idx(thisc);

[indxx, tfx] = listdlg('PromptString', {'Select groups'}, ...
    'SelectionMode', 'multiple', 'ListString', string(cLi));
if tfx == 1
    [ax, bx] = view();
    idx = ismember(ci, indxx);

    %         answer = questdlg('Select or unselect?','','Select', 'Unselect',...
    %                           'Cancel', 'Select');
    %         if strcmp(answer, 'Select')
    %               % do nothing
    %         elseif strcmp(answer, 'Unselect')
    %             idx=~idx;
    %         else
    %             return;
    %         end

    fw = gui.gui_waitbar;
    scex = selectcells(sce, idx);
    % scex.c=cLi(ci(idx));
    scex.c = ci(idx);
    gui.gui_waitbar(fw);

    hFigure = figure();
    %set(f,'WindowStyle','modal');
    %set(hFigure, 'MenuBar', 'none');
    %set(hFigure, 'ToolBar', 'none');


    UitoolbarHandle = uitoolbar('Parent', hFigure);
    % set(UitoolbarHandle, 'Tag', 'FigureToolBar', ...
    %     'HandleVisibility', 'off', 'Visible', 'on');

    pt = uipushtool(UitoolbarHandle, 'Separator', 'off');
    [img, map] = imread(fullfile(mfolder, '..', 'resources', 'list.gif'));
    ptImage = ind2rgb(img, map);
    pt.CData = ptImage;
    pt.Tooltip = 'Cluster Cells';
    pt.ClickedCallback = @clustercells;

    pt = uipushtool(UitoolbarHandle, 'Separator', 'off');
    [img, map] = imread(fullfile(mfolder, '..', 'resources', 'list2.gif'));
    ptImage = ind2rgb(img, map);
    pt.CData = ptImage;
    pt.Tooltip = 'Close Window and Return';
    pt.ClickedCallback = @returnbackparent;

    h = gui.i_gscatter3(scex.s, scex.c);
    view(ax, bx);

    waitfor(hFigure);

end

    function clustercells(~, ~)
        k = gui.i_inputnumk(5);
        try
            if ~isempty(k)
                fw2 = gui.gui_waitbar;
                newclassidenty = run.mt_SC3(scex.X, k);
                gui.gui_waitbar(fw2);
            end
            delete(h);
            h = gui.i_gscatter3(scex.s, newclassidenty);
            view(ax, bx);
        catch ME
            gui.gui_waitbar(fw2);
            errordlg(ME.message);
        end
end

        function returnbackparent(~, ~)
            if ~isempty(newclassidenty)
                requirerefresh = true;
            end
            highlightindex = idx;
            close(hFigure);
    end

    end