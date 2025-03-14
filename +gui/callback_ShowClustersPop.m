function callback_ShowClustersPop(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);

answer = gui.myQuestdlg(FigureHandle, ['Select a grouping variable and ' ...
    'show cell groups in new figures individually?']);
if ~strcmp(answer, 'Yes'), return; end

[thisc, ~] = gui.i_select1class(sce, true,'','',FigureHandle);
if isempty(thisc), return; end
[c, cL] = grp2idx(thisc);
% [c, cL, noanswer] = gui.i_reordergroups(thisc);
% if noanswer, return; end
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

    %cLa=getappdata(FigureHandle,'cL');
    %if ~isempty(cLa) && length(cL)==length(cLa)
    %    cL=cLa;
    %end
    cmv = 1:max(c);
    idxx = cmv;
    [cmx] = countmember(cmv, c);



%answer = gui.myQuestdlg(FigureHandle, 'Sort by size of cell groups?');
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

hx = gui.myFigure(FigureHandle);

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
hx.addCustomButton('off', @in_scgeatoolsce, "icon-mat-touch-app-10.gif", 'Extract and Work on Separate SCEs...');
hx.show(FigureHandle);
catch ME
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
end

    function in_scgeatoolsce(~,~)
        %answer1 = gui.myQuestdlg(FigureHandle, 'Extract cells from different groups and view new SCEs, or save new SCEs?','',...
        %    'View SCEs','Save SCEs','Cancel','View SCEs');
        answer1 = gui.myQuestdlg(FigureHandle, 'Extract cells and make new SCEs?','');
        switch answer1
            case {'Cancel','No'}
                return;
            case {'Yes','View SCEs'}
                [idx] = in_selectcellgrps(cL(idxx), FigureHandle);
                if isempty(idx), return; end
                cL2=cL(idxx);
                % currentColormap = colormap;
                % figure(FigureHandle)
                % colormap(currentColormap);
                
                s=0;
                for ik=1:length(idx)
                    % scev = SCEV{idx(ik)};
                    scev = sce.selectcells(SCEV{idx(ik)});
                    p = scgeatool(scev,'useuifig', ...
                        gui.i_isuifig(FigureHandle));
                    p.Name=matlab.lang.makeValidName(cL2{idx(ik)});
                    % p.Position([2])=p.Position([2])-s*30;
                    % p.Position([1])=p.Position([1])+s*30;
                    % p.Position([3 4])=p.Position([3 4])*0.8;
                    s=s+1;
                    pause(0.5);
                end
                if isvalid(hx)
                    hx.closeFigure;
                end
           case 'Save SCEs'
                answer2=gui.myQuestdlg(FigureHandle, 'Where to save files?','',{'Use Temporary Folder', ...
                    'Select a Folder','Cancel'},'Use Temporary Folder');
                switch answer2
                    case 'Select a Folder'
                        [seltpath] = uigetdir(deflt);
                        if isvalid(FigureHandle) && isa(FigureHandle, 'matlab.ui.Figure')
                            figure(FigureHandle);
                        end
                        
                        if seltpath==0, return; end
                        if ~isfolder(seltpath), return; end
                    case 'Use Temporary Folder'
                        seltpath = tempdir;
                    case 'Cancel'
                        return;
                end
                disp(['User selected: ', seltpath]);
                if ~isfolder(seltpath)
                    gui.myErrordlg(FigureHandle, 'Not a folder.');
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
                        answerx=gui.myQuestdlg(FigureHandle, q,'');
                    else
                        q=sprintf('Overwrite file %s?',outmatfile);
                        answerx=gui.myQuestdlg(FigureHandle, q,'');
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

function [idx] = in_selectcellgrps(grpv, FigureHandle)
    idx=[];

       if gui.i_isuifig(FigureHandle)
            [indx2, tf2] = gui.myListdlg(FigureHandle, grpv, ...
                'Select Group(s):');
        else
            [indx2, tf2] = listdlg('PromptString', ...
                {'Select Group(s):'}, ...
                'SelectionMode', 'multiple', 'ListString', grpv, ...
                'InitialValue', 1, 'ListSize', [220, 300]);
       end


    if tf2 == 1
        idx = indx2;
    end
end
