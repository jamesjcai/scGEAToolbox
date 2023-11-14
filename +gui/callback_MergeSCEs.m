function [requirerefresh, s] = callback_MergeSCEs(src, sourcetag)
requirerefresh = false;
s = "";
answer = questdlg('Current SCE will be replaced. Continue?');
if ~strcmp(answer, 'Yes'), return; end
FigureHandle = src.Parent.Parent;

keepbatchid=true;
answer = questdlg('Keep original batch IDs of cells in the input SCEs?','');
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
            warndlg('No SCE variables in Workspace.');
            return;
        elseif sum(valididx) < 2
            warndlg('Need at least two SCEs in Workspace.');
            return;
        end

        b = b(:, valididx);
        a = a(valididx);


        %sce=guidata(FigureHandle);

        [indx, tf] = listdlg('PromptString', {'Select SCEs:'}, ...
            'liststring', b(1, :), 'SelectionMode', 'multiple');
        if tf == 1
            if length(indx) < 2
                warndlg('Need at least two selected SCEs.');
                return;
            end

            answer = questdlg('Which set operation method to merge genes?', 'Merging method', ...
                'Intersect', 'Union', 'Intersect');
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
                fw = gui.gui_waitbar;
                sce = sc_mergesces(insce, methodtag, keepbatchid);
                guidata(FigureHandle, sce);
                requirerefresh = true;
            catch ME
                gui.gui_waitbar(fw, true);
                errordlg(ME.message);
                return;
            end
            gui.gui_waitbar(fw);
        else
            return;
        end
    case 2
        % warndlg("This function is under development.");

        [fname, pathname] = uigetfile({'*.mat', 'SCE Data Files (*.mat)'; ...
            '*.*', 'All Files (*.*)'}, ...
            'Select SCE Data Files', 'MultiSelect', 'on');
        if isequal(fname, 0), return; end
        if ~iscell(fname)
            errordlg("This function needs at least two SCE data files.");
                return;
            end
            answer = questdlg('Which set operation method to merge genes?', 'Merging method', ...
                'Intersect', 'Union', 'Intersect');
            if ~ismember(answer, {'Union', 'Intersect'}), return; end
            methodtag = lower(answer);

            fw = gui.gui_waitbar;
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
                gui.gui_waitbar(fw, true);
                errordlg(ME.message);
                return;
            end
            gui.gui_waitbar(fw);
end
end % end of function
