function [requirerefresh, s] = callback_MergeSCEs(src, sourcetag)
requirerefresh = false;
s = "";
[FigureHandle] = gui.gui_getfigsce(src);

answer = gui.myQuestdlg(FigureHandle, 'Current SCE will be replaced. Continue?');
if ~strcmp(answer, 'Yes'), return; end


keepbatchid=true;
answer = gui.myQuestdlg(FigureHandle, 'Keep original batch IDs of cells in the input SCEs?','');
switch answer
    case 'Yes'
        keepbatchid=true;
    case 'No'
        keepbatchid=false;
    case 'Cancel'
        return;
end


switch sourcetag
    case 1
        a = evalin('base', 'whos');
        b = struct2cell(a);
        valididx = ismember(b(4, :), 'SingleCellExperiment');
        if sum(valididx) < 1
            gui.myWarndlg(FigureHandle, 'No SCE variables in Workspace.');
            return;
        elseif sum(valididx) < 2
            gui.myWarndlg(FigureHandle, 'Need at least two SCEs in Workspace.');
            return;
        end

        b = b(:, valididx);
        a = a(valididx);


        %sce=guidata(FigureHandle);

        if gui.i_isuifig(FigureHandle)
            [indx, tf] = gui.myListdlg(FigureHandle, b(1, :), 'Select SCEs:');
        else
            [indx, tf] = listdlg('PromptString', {'Select SCEs:'}, ...
                'liststring', b(1, :), ...
                'SelectionMode', 'multiple', ...
                'ListSize', [220, 300]);
        end


        if tf == 1
            if length(indx) < 2
                gui.myWarndlg(FigureHandle, 'Need at least two selected SCEs.');
                return;
            end

            answer = gui.myQuestdlg(FigureHandle, 'Which set operation method to merge genes?', 'Merging method', ...
                {'Intersect', 'Union'}, 'Intersect');
            if ~ismember(answer, {'Union', 'Intersect'}), return; end
            methodtag = lower(answer);
            try
                insce = cell(1, length(indx));
                s = "";
                for k = 1:length(indx)
                    insce{k} = evalin('base', a(indx(k)).name);
                    s = sprintf('%s,%s', s, a(indx(k)).name);
                end
                s = s(2:end);
                fprintf('>> sce=sc_mergesces({%s},''%s'',true);\n', s, methodtag);
                fw = gui.myWaitbar(FigureHandle);
                sce = sc_mergesces(insce, methodtag, keepbatchid);
                guidata(FigureHandle, sce);
                requirerefresh = true;
            catch ME
                gui.myWaitbar(FigureHandle, fw, true);
                errordlg(ME.message);
                return;
            end
            gui.myWaitbar(FigureHandle, fw);
        else
            return;
        end
    case 2
        % gui.myWarndlg(FigureHandle, "This function is under development.");

        [fname, pathname] = uigetfile({'*.mat', 'SCE Data Files (*.mat)'; ...
            '*.*', 'All Files (*.*)'}, ...
            'Select SCE Data Files', ...
            'MultiSelect', 'on');
        figure(FigureHandle);
        if isequal(fname, 0), return; end
        if ~iscell(fname)
            errordlg("This function needs at least two SCE data files.");
            return;
        end

        answer = gui.myQuestdlg(FigureHandle, 'Which set operation method to merge genes?', 'Merging method', ...
            {'Intersect', 'Union'}, 'Intersect');
        if ~ismember(answer, {'Union', 'Intersect'}), return; end
        methodtag = lower(answer);

        fw = gui.myWaitbar(FigureHandle);
        try
            scelist = cell(length(fname));
            s = "";
            for k = 1:length(fname)
                scefile = fullfile(pathname, fname{k});
                load(scefile, 'sce');
                sce.metadata = [sce.metadata; fname{k}];
                scelist{k} = sce;
                s = sprintf('%s, %s', s, fname{k});
            end
            s = s(2:end);
            pause(1)
            sce = sc_mergesces(scelist, methodtag, keepbatchid);
            guidata(FigureHandle, sce);
            requirerefresh = true;
        catch ME
            gui.myWaitbar(FigureHandle, fw, true);
            errordlg(ME.message);
            return;
        end
        gui.myWaitbar(FigureHandle, fw);
end
end % end of function
