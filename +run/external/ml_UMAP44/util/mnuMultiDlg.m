function [idxs, cancelled, mainDlg, raIfRemembered, checkBoxes, allChb, sortGui]=...
    mnuMultiDlg(args, title, options, defaults, singleOnly, ...
    radioOrCheckBox, xtraCmp1, xtraWhere1, xtraCmp2, xtraWhere2, ...
    itemsPerScrollWindow,  where, subTitle)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

%   Since 2013 this module has evolved like a platypus.  Confusing.
%   Somewhat simplified by Gui.Ask but truly deserves a new wrapper
%   with addParameter() argument preparation.  The module is directly
%   called too many times to put in the trash however.  Sorry.

cancelled=true;
raIfRemembered=[];
allChb=[];
checkBoxParent=[];
checkBoxParentCount=0;
if nargin<13
    subTitle=[];
    if nargin<12
        where='center';
        if nargin<11
            itemsPerScrollWindow=11;
            if nargin <10
                xtraWhere2='South';
                if nargin<9
                    xtraCmp2=[];
                    if nargin<8
                        xtraWhere1='South';
                        if nargin<7
                            xtraCmp1=[];
                            if nargin<6
                                radioOrCheckBox=false;
                                if nargin<5
                                    singleOnly=false;
                                    if nargin<4
                                        defaults=[];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
nOptions=length(options);
checkBoxes=[];
noCancel=isstruct(args) && isfield(args, 'noCancel') && args.noCancel;
disableChoices=isstruct(args) && isfield(args, 'disableChoices') ...
    && args.disableChoices;
if isstruct(args) && isfield(args, 'disableIdxs') 
    disableIdxs=args.disableIdxs;
else
    disableIdxs=[];
end
if isstruct(args) && isfield(args, 'allMsg')
    allMsg=args.allMsg;
elseif nOptions>1
    allMsg='All';
end
[guide, where, property, properties, defaults, icon, javaWindow, ...
    choiceTitle, checkFnc, modal, ~, ~, rememberId, ~, ...
    checkBoxFnc, items]=...
    decodeMsg(args, defaults, where);
if items>2
    itemsPerScrollWindow=items;
end
if ~isempty(rememberId)
    ra=RememberedAnswers;
    idxs=ra.getIdxs(rememberId);
    if ra.isRemembered(rememberId)
        mainDlg=[];
        cancelled=false;
        raIfRemembered=ra;
        return;
    end
end
cancelled=true;
if ischar(options{1})
    jsa=javaArray('java.lang.String', nOptions);
    for ix=1:nOptions
        jsa(ix)=java.lang.String(options{ix});
    end
else
    jsa=javaArray('java.lang.Object', nOptions);
    for iy=1:nOptions
        jsa(iy)=options{iy};
    end
end
hasIcon=isempty(icon) || ~strcmp('none', icon);
pnl=Gui.BorderPanel(1,1);
if hasIcon
        pnl.setBorder(javax.swing.BorderFactory.createEmptyBorder(4, 7, 4, 5));
    if isempty(icon)
        pnl.add(javaObjectEDT('javax.swing.JLabel', Gui.Icon('facs.gif')), 'West');
    else
        pnl.add(javaObjectEDT('javax.swing.JLabel', Gui.Icon(icon)), 'West');
    end
    pnl.add(javaObjectEDT('javax.swing.JLabel', Gui.Icon('blank.png')), 'East');
else
    pnl.setBorder(javax.swing.BorderFactory.createEmptyBorder(2, 4, 2, 4));
end
if jsa.length<itemsPerScrollWindow
    itemsPerScrollWindow=jsa.length;
end
if itemsPerScrollWindow>18
    itemsPerScrollWindow=18;
end
if radioOrCheckBox
    dfltIdxs=0;
    if ~isempty(defaults)
        if singleOnly
            dfltIdxs=defaults(1);
        else
            dfltIdxs=defaults+1;
        end
    end
    if isfield(args, 'sortProps') 
        [jsc,bg,checkBoxes, innerPnl]=Radio.Panel4(jsa, dfltIdxs, ...
            itemsPerScrollWindow, ~singleOnly, false);
    else
        radioWithExtraComponents=isstruct(args)...
            && isfield(args, 'radioWithExtraComponents') ...
            && args.radioWithExtraComponents;
        [jsc,bg,checkBoxes, innerPnl]=Radio.Panel4(jsa, dfltIdxs, ...
            itemsPerScrollWindow, ~singleOnly, true, ...
            radioWithExtraComponents);
    end
    if ~singleOnly && isa(jsa(1),'javax.swing.JPanel')
        priorCmps=Gui.GetJavaComponents(innerPnl);
    end
    if ~isempty(choiceTitle)
        Gui.SetTitledBorder(choiceTitle, jsc);
    end
    innerPnl.setBackground(java.awt.Color(1, 1, 1))
    
    if ~isempty(disableIdxs)
        try
            nDisableIdxs=length(disableIdxs);
            for iii=1:nDisableIdxs
               di=disableIdxs(iii);
               checkBoxes.get(di-1).setEnabled(false);
            end
        catch ex
            ex.getReport
        end
    elseif disableChoices
        nCheckBoxes=checkBoxes.size;
        for iii=0:nCheckBoxes-1
            checkBoxes.get(iii).setEnabled(false);
        end
    end
    maxSelections=0; minSelections=1;
    if ~singleOnly && isstruct(args) 
        initCheckBoxes;
        if isfield(args, 'max')
            maxSelections=args.max;
            MatBasics.RunLater(@(h,e)checkSelections(), 1);
        end
        if isfield(args, 'min')
            minSelections=args.min;
            MatBasics.RunLater(@(h,e)checkSelections(), 1);
        end
    end
else
    lst=javaObjectEDT('javax.swing.JList', jsa);
    if disableChoices
        lst.setEnabled(false);
    end
    lst.setBorder(javax.swing.BorderFactory.createEmptyBorder(4, 10, 4, 8));
    jListbox = handle(lst, 'CallbackProperties');    
    set(jListbox, 'MousePressedCallback',@myCallbackFcn);
    % Define the mouse-click callback function
    doubleClicked=false;
    scroll=javaObjectEDT('javax.swing.JScrollPane', lst);
    limit=itemsPerScrollWindow;
    lst.setVisibleRowCount(limit);
    d=lst.getPreferredSize;
    if d.width>650
        d.width=650;
        lst.setPreferredSize(d);
    end
    jsc=scroll;
end
if ~isempty(xtraCmp1) && ~strcmp(xtraWhere1,'south buttons')...
        && ~strcmp(xtraWhere1,'south west buttons')
    pnl2=Gui.BorderPanel(2,8);
    pnl2.add(xtraCmp1, xtraWhere1);
    if ~isempty(xtraCmp2)
        pnl2.add(xtraCmp2, xtraWhere2);
    end
    pnl.add(pnl2, 'Center');
    pnl3=Gui.BorderPanel(2,10);
    pnl2.add(pnl3, 'Center');
    pnl2=pnl3;
else
    pnl2=pnl;
end
allChbPnl=[];
flowPnl=Gui.FlowPanelCenter;
flowPnl.add(jsc);
if ~singleOnly && radioOrCheckBox
    pnl3=Gui.BorderPanel(0,0);
    noAll=isstruct(args) && isfield(args, 'noAll') && args.noAll;
    if ~noAll && nOptions>1
        allChb=javaObjectEDT('javax.swing.JCheckBox', allMsg);
        allChb.setMnemonic('a');
        if length(unique(defaults))==length(options)
            allChb.setSelected(true);
        end
        allH=handle(allChb,'CallbackProperties');
        allChbPnl=Gui.BorderPanel(0,0);
        allChbPnl.add(allChb, 'West');
        pnl3.add(allChbPnl, 'North');
        set(allH, 'ActionPerformedCallback', @(h,e)doAll(h,e));
    end
    pnl3.add(flowPnl, 'Center');
    checkBoxPnl=pnl3;
    flowPnl=Gui.FlowPanelCenter;
    flowPnl.add(pnl3);
    pnl2.add(flowPnl, 'Center');
else
    pnl2.add(flowPnl, 'Center');
end
pnl3=Gui.BorderPanel(0,0);
if hasIcon
    pnl3.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 28, 0, 0));
end
if ischar(guide)
    lbl=javaObjectEDT('javax.swing.JLabel', guide);
    lbl.setHorizontalAlignment(javax.swing.JLabel.CENTER);
    pnl3.add(lbl, 'North');
else
    pnl3.add(guide, 'North');
end
pnl2.add(pnl3, 'North');
if ~isempty(subTitle) && ischar(subTitle)
    lbl2=javaObjectEDT('javax.swing.JLabel', subTitle);
    lbl2.setHorizontalAlignment(javax.swing.JLabel.CENTER);
    pnl.add(lbl2, 'North');
end
if nargin<=3 || isempty(defaults)
    defaults=0;
end
if ~radioOrCheckBox
    if singleOnly
        lst.setSelectionMode(1);
    end
    if ~isempty(defaults)
        lst.setSelectedIndices(int32(defaults));
        lst.ensureIndexIsVisible(defaults(1));
    end
end

if isempty(javaWindow)
    jFrame = Gui.ParentFrame;
else
    jFrame = javaWindow;
end
mainDlg=javaObjectEDT('javax.swing.JDialog', jFrame);
mainDlg.setBackground(java.awt.Color.WHITE);
dlg=handle(mainDlg, 'CallbackProperties');
set(dlg, 'WindowClosingCallback', @(h,e)windowClose());
if ~isempty(title)
    mainDlg.setTitle(title);
end
if ~noCancel
    done=javaObjectEDT('javax.swing.JButton', 'Ok');
else
    done=javaObjectEDT('javax.swing.JButton', 'Done');
end
doneH=handle(done,'CallbackProperties');
set(doneH, 'ActionPerformedCallback', @(h,e)close(true));
cancel=javaObjectEDT('javax.swing.JButton', 'Cancel');
cancel.setIcon(Gui.Icon('cancel.gif')); 
cancelH=handle(cancel,'CallbackProperties');
set(cancelH, 'ActionPerformedCallback', @(h,e)close(false));
try
    edu.stanford.facs.swing.CpuInfo.registerEscape(mainDlg, cancel);
catch
end
south=Gui.BorderPanel(0,0);
southSouth=Gui.BorderPanel(0,0);
southEast=Gui.FlowPanel;
southSouth.add(southEast, 'South');
if ~isempty(rememberId)
    rememberCb=RememberedAnswers.GetCheckBox;
    southEast.add(rememberCb);
    southEast.add(Gui.Label(' '));
end
if ~noCancel
    southEast.add(cancel);
end
southEast.add(done);
south.add(southSouth, 'East');
if ~isempty(xtraCmp1) && strcmp(xtraWhere1,'south west buttons') 
    south.add(xtraCmp1, 'West');  
    if ~isempty(xtraCmp2)
        if strcmp('South', xtraWhere2)
            south.add(xtraCmp2, 'North');
        else
            pnl.add(xtraCmp2, xtraWhere2);
        end
    end
else
    if ~isempty(xtraCmp1) && strcmp(xtraWhere1,'south buttons')
        south.add(xtraCmp1, 'Center');
    end
end
mainDlg.getRootPane.setDefaultButton(done);
pnl.add(south, 'South');
pnl.setBorder(javax.swing.BorderFactory.createEmptyBorder(2, 5, 10, 5));
mainDlg.add(pnl);
sortGui=[];
if radioOrCheckBox 
    if isempty(allChbPnl)
        allChbPnl=Gui.BorderPanel(0,0);
        pnl3.add(allChbPnl, 'South');
        allMsg='All';
        allChb=[];
    end
    if (isfield(args, 'sortGuiAlways') && ...
            args.sortGuiAlways) ...
            || (singleOnly && nOptions>7) ...
            || (~singleOnly && nOptions>2)
        sortGui=SortGui(mainDlg, allChb, allMsg, allChbPnl, options, ...
            checkBoxes, innerPnl);
        if isa(jsa(1),'javax.swing.JPanel')
            sortGui.fncRefresh=@refreshPanelOrder;
        end
        if isfield(args, 'sortProps') 
            if isfield(args, SortGui.PROP_DEFAULT_IDX)
                searchDfltIdx=args.sortDefaultIdx;
                searchOn=true;
            else
                searchDfltIdx=0;
                searchOn=isfield(args, SortGui.PROP_SEARCH2) ...
                    && args.sortSearch;
            end
            sortGui.setProperties(args.sortProps, ...
                args.sortProp, searchOn, searchDfltIdx);
        end
    end
else
    sortGui=[];
end
mainDlg.pack;
if ~isempty(jFrame)
    mainDlg.setLocationRelativeTo(jFrame);
end
Gui.LocateJava(mainDlg, javaWindow, where);
if radioOrCheckBox
    nCh=checkBoxes.size;
    for iCh=1:nCh
        cb1=checkBoxes.get(iCh-1);
        cb2= handle(cb1, 'CallbackProperties');
        set(cb2, 'ActionPerformedCallback', @(h,e)innerChbCb(h,e));
    end
end
cancelled=true;
if ~ispc
    setAlwaysOnTopTimer(mainDlg)
end
MatBasics.RunLater(@(h,e)noTip(), .25)
if ~isempty(sortGui)
    sortGui.setAllChbText;
end
lastPick=[];
mainDlg.setModal(modal);
Gui.SetJavaVisible(mainDlg);
drawnow;
conclude;

    function noTip
        BasicMap.Global.closeToolTip;
    end

    function conclude
        if ~radioOrCheckBox
            javaIdxs=lst.getSelectedIndices;
            N=length(javaIdxs);
        else
            doubleClicked=false;
            N=0;
        end
        idxs=[];
        if doubleClicked || ~cancelled
            if N==0
                if radioOrCheckBox
                    if singleOnly
                        idxs=Radio.Choice(bg);
                    else
                        idxs=getSelectedIdxs;
                    end
                else
                    %idxs=1;
                end
            else
                for i=1:N
                    idx=javaIdxs(i)+1;
                    idxs(end+1)=idx;
                    answer=StringArray.IndexOf(options, idx);
                    disp(answer);
                end
            end
        end
        if ~cancelled
            if ~isempty(property) && ~isempty(properties)
                if singleOnly
                    properties.set(property, num2str(idxs));
                else
                    properties.set(property, num2str(idxs-1));
                end
            end
            if ~isempty(rememberId)
                if rememberCb.isSelected
                    rememberAnswer(idxs);
                end
            end
        end
    end

    function windowClose
        close(false);
    end

    function doAll(h,e)
        isSelected=h.isSelected;
        N2=checkBoxes.size;
        for ii=1:N2
            cb=checkBoxes.get(ii-1);
            if ~isempty(Gui.WindowAncestor(cb))
                cb.setSelected(isSelected);
            end
        end
        innerChbCb(h,e);
    end

    function close(saved)
        if isstruct(args) && isfield(args, 'closeFnc')
            if ~isempty(args.closeFnc)
                feval(args.closeFnc, saved, idxs, checkBoxes);
                return;
            end
        end
        cancelled=~saved;
        if ~isempty(checkFnc)
            conclude;
            ok=feval(checkFnc, idxs, cancelled, mainDlg);
            if ~ok
                return;
            end
        end
        mainDlg.dispose;
    end

    function myCallbackFcn(~,jEventData)
        % Determine the click type
        % (can similarly test for CTRL/ALT/SHIFT-click)
        if jEventData.getClickCount==2
            w=Gui.Wnd(lst);
            if ~isempty(w)
                doubleClicked=true;
                if modal
                    w.dispose;
                else
                    close(true);
                end
            end
        end
    end

    function rememberAnswer(idxs)
        curFig=get(0, 'currentFigure');
        [~, ~, ~, ~, quadrant]=Gui.FindScreen(curFig);
        ra=RememberedAnswers;
        if strcmpi('west', quadrant{2})
            where='north east++';
        else
            where='north west++';
        end
        Gui.LocateJava(mainDlg, javaWindow, where);
        mainDlg.setModal(false);
        mainDlg.setVisible(true);
        
        ra.remember(rememberId, rememberCb, idxs);
        raIfRemembered=ra;
    end

    function [idxs_, N]=getSelectedIdxs
        [idxs_, N]=Gui.GetSelectedChbIdxs(checkBoxes);
    end

    function checkSelections(idxs_)
        if nargin<1
            if exist('sortGui', 'var') && ~isempty(sortGui)
                idxs_=sortGui.setAllChbText;
            else
                idxs_=Gui.GetSelectedChbIdxs(checkBoxes);
            end
        end
        if maxSelections>0 || minSelections>0
            nIdxs=length(idxs_);
            if (maxSelections>0 && nIdxs>maxSelections) ...
                || nIdxs < minSelections
                mx=maxSelections;
                if maxSelections>0 && nIdxs>maxSelections
                    mxSuffix=sprintf(['<li><i>Remove %d of ' ...
                        'your %d selections....</i>'], ...
                        nIdxs-maxSelections, nIdxs);
                else
                    if maxSelections==0
                        mx=checkBoxes.size;
                    end
                    mxSuffix='';
                end
                if nIdxs<minSelections
                    mnSuffix=sprintf('<li>Select at least %s', ...
                        String.Pluralize2('item', minSelections));
                else
                    mnSuffix='';
                end
                if nIdxs>0
                    cmp=checkBoxes.get(idxs_(end)-1);
                else
                    cmp=innerPnl;
                end
                if minSelections==maxSelections
                    prefix='exactly';
                else
                    prefix=sprintf('between %d and', minSelections);
                end
                tip=Html.Sprintf([' %s&nbsp;&nbsp;Pick %s %d items!' ...
                    '<ul>%s%s</ul>'], ...
                    Html.ImgXy('error.png', [], ...
                    BasicMap.Global.adjustHighDef(.4)),...
                    prefix, mx, mnSuffix, mxSuffix);
                innerPnl.setBorder(...
                    javax.swing.BorderFactory.createLineBorder...
                    (java.awt.Color.RED, 4));
                Gui.Shake(cmp, 4, tip);
            else
                innerPnl.setBorder([]);
            end
        end
    end

    function initCheckBoxes
        try
            checkBoxParent=checkBoxes.get(0).getParent;
            checkBoxParentCount=checkBoxParent.getComponentCount;
        catch
        end
    end
    
    function selectGroup(curPick)
        if isempty(lastPick)
            return;
        end
        if checkBoxParent==0
            initCheckBoxes;
        end
        if ~lastPick.isSelected
            return;
        end
        for ii=0:checkBoxParentCount-1
            cb=checkBoxParent.getComponent(ii);
            if isequal(cb, lastPick)
                for iix=ii+1:checkBoxParentCount-1
                    cb=checkBoxParent.getComponent(iix);
                    if isequal(cb, curPick)
                        break;
                    end
                    cb.setSelected(true);
                end
                break;
            elseif isequal(cb, curPick)
                for iix=ii+1:checkBoxParentCount-1
                    cb=checkBoxParent.getComponent(iix);
                    if isequal(cb, lastPick)
                        break;
                    end
                    cb.setSelected(true);
                end
                break;
            end
        end
    end

    function innerChbCb(h, e)
        modifiers=e.getModifiers;
        if ~isempty(sortGui)
            idxs_=sortGui.setAllChbText;
        else
            idxs_=Gui.GetSelectedChbIdxs(checkBoxes);
        end
        if ~singleOnly
            cb=e.getSource;
            if ~isequal(cb, allChb)
                if cb.isSelected
                    if modifiers==17 || modifiers==21
                        selectGroup(cb);
                        sortGui.setAllChbText
                    end
                    lastPick=cb;
                end
            else
                lastPick=[];
            end
            checkSelections(idxs_);
        end
        if ~isempty(checkBoxFnc)
            feval(checkBoxFnc, h, e, idxs_, checkBoxes);
        end
    end

    function innerPnl=refreshPanelOrder
        checkBoxPnl.remove(jsc);
        [jsc,~,~, innerPnl]=Radio.Panel2(jsa, [], ...
            itemsPerScrollWindow, true, [], priorCmps, ...
            sortGui.sortIdxs, sortGui.visibleIdxs);
        checkBoxPnl.add(jsc, 'Center');
        jsc.repaint;
    end

end
