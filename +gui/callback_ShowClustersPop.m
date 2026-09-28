function callback_ShowClustersPop(src, ~)


[FigureHandle, sce] = gui.gui_getfigsce(src);
% The app, kept apart from SRC: IN_CALLBACK_SCGEATOOLSCE takes a SRC of its
% own, and a nested function's argument is this workspace's variable, so the
% first button click would overwrite it with the button.
appsrc = src;
if ~isempty(FigureHandle) && pkg.i_isvalid(FigureHandle) && FigureHandle.Visible == "on"
    figure(FigureHandle);
    cleanupObj = onCleanup(@() gui.i_raisefig(FigureHandle));
end

answer = gui.myQuestdlg(FigureHandle, ['Select a grouping variable and ' ...
'show cell groups in new figures individually?']);
if ~strcmp(answer, 'Yes'), return; end

% Several grouping variables may be picked; they cross into one composite
% label per cell ("Macrophages | IL"). Downstream treats thisc as a
% per-cell label vector, so the composite needs no special handling.
[thisc, ~] = gui.i_selectnclass(sce, true,'','',FigureHandle);
if isempty(thisc), return; end
[c, cL] = findgroups(string(thisc));
if max(c)==1
    gui.myHelpdlg(FigureHandle, sprintf('Only one type of cells: %s',cL{1}))
    return;
end


   SCEV=cell(max(c),1);
   try
        for k=1:max(c)
            SCEV{k} = c==k;
        end
    catch ME
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end

cmv = 1:max(c);
idxx = cmv;
[cmx] = countmember(cmv, c);


[~, idxx] = sort(cmx, 'descend');
SCEV = SCEV(idxx);

try
    sces = sce.s;
    % The app's own scatter: findall returned every scatter on the window,
    % and h.ZData errored for none or for more than one.
    if isa(src, 'matlab.apps.AppBase') && pkg.i_isvalid(src.h)
        h = src.h;
    else
        h = findall(FigureHandle, 'type', 'scatter');
    end
    if isempty(h) || isempty(h(1).ZData), sces = sce.s(:, 1:2); end

    [para] = gui.i_getoldsettings(src, FigureHandle);

    totaln = max(c);
    numfig = ceil(totaln/9);

    % -------------

    hx = gui.myFigure(FigureHandle, true);

    tabgp = uitabgroup(hx.FigHandle);
    for nf = 1:numfig
        tab{nf} = uitab(tabgp, 'Title', sprintf('Tab%d',nf));
        axes('parent',tab{nf});
        for k=1:9
            kk = (nf - 1) * 9 + k;
            if kk <= totaln
                ax{nf, k} = subplot(3,3,k);
                gui.i_gscatter3(sces, c, 3, cmv(idxx(kk)));
                set(ax{nf, k}, 'XTick', []);
                set(ax{nf, k}, 'YTick', []);
                b = cL{idxx(kk)};
                title(ax{nf, k}, strrep(b, '_', "\_"));
                a = sprintf('%s (%.2f%%)', ...
                    pkg.i_plural(cmx(idxx(kk)), 'cell'), ...
                    100*cmx(idxx(kk))/length(c));
                fprintf('%s in %s\n', a, b);
                subtitle(ax{nf, k}, a);
                box(ax{nf, k}, 'on');
            end
        end
        % The figure, so every panel gets it: on hx.AxHandle, the hidden
        % axes under the tab group, it reached none of them.
        if isfield(para, 'oldColorMap'), colormap(hx.FigHandle, para.oldColorMap); end
    end
    hx.addCustomButton('off', @in_callback_scgeatoolsce, "icon-mat-touch-app-10.gif", 'Work on a Group, or Save Groups as SCEs...');
    hx.show(FigureHandle);
catch ME
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
end

function in_callback_scgeatoolsce(~, ~)
        figure(hx.FigHandle);

        % Working on a group replaces the dataset in the main window, with
        % an Undo snapshot, rather than opening a scgeatoolApp window per
        % group over it. Keeping several groups as separate datasets is what
        % Save SCEs is for. (This dialog used to offer only Yes/No/Cancel,
        % so Save SCEs could not be reached.)
        optHere = 'Work on One Group';
        optSave = 'Save SCEs';
        answer1 = gui.myQuestdlg(hx.FigHandle, ['Work on one group in the ', ...
            'main window (Edit > Undo brings the rest back), or save groups ', ...
            'as separate SCE files?'], '', {optHere, optSave, 'Cancel'}, optHere);
        switch answer1
            case optHere
                [idx] = in_selectcellgrps(cL(idxx), hx.FigHandle, false);
                if isempty(idx), return; end
                scev = copy(sce).selectcells(SCEV{idx}); % OK
                if isa(appsrc, 'matlab.apps.AppBase')
                    % Before the data changes, not in the menu handler:
                    % only this button changes it, and a snapshot for the
                    % display alone would push out the last real Undo.
                    gui.i_snapshot(appsrc, 'Work on Cell Group');
                    hx.closeFigure;
                    gui.i_replacesce(appsrc, scev, sce.NumCells);
                else
                    scgeatool(scev);
                end
            case optSave
                answer2=gui.myQuestdlg(hx.FigHandle, 'Where to save files?','',{'Use Temporary Folder', ...
                    'Select a Folder','Cancel'},'Use Temporary Folder');
                switch answer2
                    case 'Select a Folder'
                        if gui.i_isuifig(hx.FigHandle)
                            [seltpath] = uigetdir(hx.FigHandle);
                        else
                            [seltpath] = uigetdir();
                        end

                        if seltpath==0, return; end
                        if ~isfolder(seltpath), return; end
                    case 'Use Temporary Folder'
                        seltpath = tempdir;
                    otherwise
                        return;
                end
                disp(['User selected: ', seltpath]);
                if ~isfolder(seltpath)
                    gui.myErrordlg(hx.FigHandle, 'Not a folder.');
                    return;
                end

                [idx] = in_selectcellgrps(cL(idxx), hx.FigHandle, true);
                cL2=cL(idxx);
                if isempty(idx), return; end
                for ik=1:length(idx)
                    scev = copy(sce).selectcells(SCEV{idx(ik)}); % OK

                    scev=scev.qcfilter;
                    outmatfile=sprintf('%s.mat', ...
                        matlab.lang.makeValidName(cL2{idx(ik)}));
                    outmatfile=fullfile(seltpath,outmatfile);
                    if ~exist(outmatfile,"file")
                        q=sprintf('Save file %s?',outmatfile);
                    else
                        q=sprintf('Overwrite file %s?',outmatfile);
                    end
                    if ~strcmp(gui.myQuestdlg(hx.FigHandle, q,''), 'Yes'), return; end
                    in_savesce(outmatfile, scev);
                end
            otherwise
                return;
        end
    end
end

function in_savesce(outmatfile, sce)
% A local function, so the variable saved as SCE is not the caller's SCE:
% assigning SCE inside the nested callback replaced the dataset it went on
% to extract the next group from.
save(outmatfile, 'sce', '-v7.3');
end

function [idx] = in_selectcellgrps(grpv, FigureHandle, multiple)
idx=[];
if multiple
    prompt = 'Select Group(s):';
    selmode = 'multiple';
else
    prompt = 'Select a Group:';
    selmode = 'single';
end

   if gui.i_isuifig(FigureHandle)
        [indx2, tf2] = gui.myListdlg(FigureHandle, grpv, ...
            prompt, grpv(1), multiple);
    else
        [indx2, tf2] = listdlg('PromptString', ...
            {prompt}, ...
            'SelectionMode', selmode, 'ListString', grpv, ...
            'InitialValue', 1, 'ListSize', [220, 300]);
   end


if tf2 == 1
    idx = indx2;
end
end
