function callback_ShowClustersPop(src, ~)
answer = questdlg('Select a grouping variable and show cell groups in new figures individually?');
if ~strcmp(answer, 'Yes'), return; end

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
[thisc, ~] = gui.i_select1class(sce);
if isempty(thisc), return; end
[c, cL] = grp2idx(thisc);
% [c, cL, noanswer] = gui.i_reordergroups(thisc);
% if noanswer, return; end
if max(c)==1
    waitfor(helpdlg(sprintf('Only one type of cells: %s',cL{1}),''))
    return;
end

%    fw = gui.gui_waitbar_adv;
    SCEV=cell(max(c),1);

   try
        for k=1:max(c)
            % gui.gui_waitbar_adv(fw, ...
            %     (k-1)/max(c), ...
            %     sprintf('Processing %s ...', cL{k}));
            % SCEV{k}=sce.selectcells(c==k);
            SCEV{k} = c==k;
        end
    catch ME
       % gui.gui_waitbar_adv(fw);
        errordlg(ME.message);
        return;
    end

    %cLa=getappdata(FigureHandle,'cL');
    %if ~isempty(cLa) && length(cL)==length(cLa)
    %    cL=cLa;
    %end
    cmv = 1:max(c);
    idxx = cmv;
    [cmx] = countmember(cmv, c);

%gui.gui_waitbar_adv(fw);

%answer = questdlg('Sort by size of cell groups?');
%if strcmpi(answer, 'Yes')
    [~, idxx] = sort(cmx, 'descend');
    SCEV = SCEV(idxx);
%end

try
    sces = sce.s;
    h = findall(FigureHandle, 'type', 'scatter');
    if isempty(h.ZData), sces = sce.s(:, 1:2); end

    [para] = gui.i_getoldsettings(src);
    totaln = max(c);
    numfig = ceil(totaln/9);

% -------------

hFig = figure('visible', 'off','Position',FigureHandle.Position);
tabgp = uitabgroup();
for nf = 1:numfig
    tab{nf} = uitab(tabgp, 'Title', sprintf('Tab%d',nf));
    axes('parent',tab{nf});
    for k=1:9
        kk = (nf - 1) * 9 + k;
        if kk <= totaln
        ax{nf,k} = subplot(3,3,k);
        gui.i_gscatter3(sces, c, 3, cmv(idxx(kk)));
        set(ax{nf, k}, 'XTick', []);
        set(ax{nf, k}, 'YTick', []);
        b = cL{idxx(kk)};
        title(strrep(b, '_', "\_"));
        a = sprintf('%d cells (%.2f%%)', ...
            cmx(idxx(kk)), ...
            100*cmx(idxx(kk))/length(c));
        fprintf('%s in %s\n', a, b);
        subtitle(a);
        box on
        end
    end
    colormap(para.oldColorMap);
end
tb = findall(hFig, 'Tag', 'FigureToolBar'); % get the figure's toolbar handle
pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
pkg.i_addbutton2fig(tb, 'off', @in_scgeatoolsce, "icon-mat-touch-app-10.gif", 'Extract and Work on Separate SCEs...');
drawnow;
set(hFig, 'visible', 'on');

%{ 
--------------------------------------------    
    for nf = 1:numfig
        f = figure('visible', 'off');
        for k = 1:9
            kk = (nf - 1) * 9 + k;
            if kk <= totaln
                %subplot(3, 3, k);
                nexttile;
                gui.i_gscatter3(sces, c, 3, cmv(idxx(kk)));
                set(gca, 'XTick', []);
                set(gca, 'YTick', []);
                b = cL{idxx(kk)};
                title(strrep(b, '_', "\_"));
                a = sprintf('%d cells (%.2f%%)', ...
                    cmx(idxx(kk)), ...
                    100*cmx(idxx(kk))/length(c));
                fprintf('%s in %s\n', a, b);
                subtitle(a);
                box on
            end
            colormap(para.oldColorMap);
        end
        P = get(f, 'Position');
        set(f, 'Position', [P(1) - 20 * nf, P(2) - 20 * nf, P(3), P(4)]);        
        tb = uitoolbar(f);
        pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
        if nf==1
            pkg.i_addbutton2fig(tb, 'off', @in_scgeatoolsce, "icon-mat-touch-app-10.gif", 'Extract and Work on Separate SCEs...');
        end
        drawnow;
        set(f, 'visible', 'on');
    end
%} 

catch ME
    errordlg(ME.message);
end

    function in_scgeatoolsce(~,~)
        %answer1 = questdlg('Extract cells from different groups and view new SCEs, or save new SCEs?','',...
        %    'View SCEs','Save SCEs','Cancel','View SCEs');
        answer1 = questdlg('Extract cells and make new SCEs?','');
        switch answer1
            case {'Cancel','No'}
                return;
            case {'Yes','View SCEs'}
                [idx] = in_selectcellgrps(cL(idxx));                
                if isempty(idx), return; end 
                for ik=1:length(idx)
                    % scev = SCEV{idx(ik)};
                    scev = sce.selectcells(SCEV{idx(ik)});
                    scgeatool(scev);
                    pause(0.5);
                end
           case 'Save SCEs'
                answer2=questdlg('Where to save files?','','Use Temporary Folder', ...
                    'Select a Folder','Cancel','Use Temporary Folder');
                switch answer2
                    case 'Select a Folder'
                        [seltpath] = uigetdir(deflt);
                        if seltpath==0, return; end
                        if ~isfolder(seltpath), return; end
                    case 'Use Temporary Folder'
                        seltpath = tempdir;
                    case 'Cancel'
                        return;
                end
                disp(['User selected: ', seltpath]);
                if ~isfolder(seltpath)
                    errordlg('Not a folder.');
                    return;
                end
               
                [idx] = in_selectcellgrps(cL(idxx));
                cL2=cL(idxx);
                if isempty(idx), return; end 
                for ik=1:length(idx)
                    % scev=SCEV{idx(ik)};
                    scev = sce.selectcells(SCEV{idx(ik)});
                    
                    scev=scev.qcfilter;
                    outmatfile=sprintf('%s.mat', ...
                        matlab.lang.makeValidName(cL2{idx(ik)}));
                    outmatfile=fullfile(seltpath,outmatfile);
                    if ~exist(outmatfile,"file")
                        q=sprintf('Save file %s?',outmatfile);
                        answerx=gui.questdlg_timer(15,q,'','Yes','No','Cancel','Yes');
                    else
                        q=sprintf('Overwrite file %s?',outmatfile);
                        answerx=questdlg(q,'');
                    end
                    switch answerx
                        case 'Yes'
                            sce = scev;
                            save(outmatfile, 'sce', '-v7.3');
                        otherwise
                            return;
                    end
                    pause(0.5);
                end
            otherwise
                return;
        end
    end
end

function [idx] = in_selectcellgrps(grpv)
    idx=[];
    [indx2, tf2] = listdlg('PromptString', ...
    {'Select Group(s):'}, ...
    'SelectionMode', 'multiple', 'ListString', grpv, ...
    'InitialValue', 1, 'ListSize', [220, 300]);
    if tf2 == 1
        idx = indx2;
    end
end
