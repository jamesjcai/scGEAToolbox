classdef SuhTree < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties
        tree; %uitree instance
        jtree; %javax.swing.JTree instance
        props;
        fncNewChildren;
        fncGetChildren;
        fncNodeExists;
        fncGetPath;
        uiNodes;
        root;
        remembered=false;
        container;
        isCollapsing=false;
        isCtrlDown=false;
        isAltDown=false;
        isMetaDown=false;
        isShiftDown=false;
        noGroupAction=false;
    end
    
    methods
        function this=SuhTree(tree, container, root, fncGetPath, ...
                fncNodeExists, fncNewChildren, fncGetChildren)
            this.tree=tree;
            this.root=root;
            this.container=container;
            this.jtree=javaObjectEDT(handle( ... %double EDT protection?
                tree.getTree,'CallbackProperties'));
            this.fncNodeExists=fncNodeExists;
            this.fncGetPath=fncGetPath;
            this.fncNewChildren=fncNewChildren;
            this.fncGetChildren=fncGetChildren;
            this.uiNodes=java.util.TreeMap;
            this.jtree.setToggleClickCount(0);
            this.setMultipleSelection;
            this.setDefaultExpandedCallback;
            set(this.tree, 'NodeCollapsedCallback', ...
                @(h,e)nodeCollapsed(this, e));
            set(this.jtree, 'KeyPressedCallback', ...
                @(src, evd)keyPressedFcn(this, evd));
            set(this.jtree, 'KeyReleasedCallback', ...
                @(src, evd)keyReleasedFcn(this, evd));
            this.jtree.setRowHeight(0);
            if BasicMap.Global.highDef
                javaMethodEDT('setTreeRowBorder', ...
                    'edu.stanford.facs.swing.SwingUtil', this.jtree, ...
                    3, 2, 2, 2);
            else
                javaMethodEDT('setTreeRowBorder', ...
                    'edu.stanford.facs.swing.SwingUtil', this.jtree, 2, 1, 1, 1);
            end
        end
        
        function suppressGroupAction(this)
            this.noGroupAction=true;
            MatBasics.RunLater(@(h,e)off(), 2.3);
            function off
                this.noGroupAction=false;
            end
        end
        
        function groupAction=isGroupAction(this)
            if this.noGroupAction
                groupAction=false;
            else
                if ispc
                    groupAction=this.isCtrlDown;
                else
                    groupAction=this.isMetaDown;
                end
            end
        end

        function keyPressedFcn(this, evd)
            this.trackStateKeys(evd);
        end
        
        function keyReleasedFcn(this, evd)
            this.trackStateKeys(evd);
        end
        
        function trackStateKeys(this, evd)
            this.isCtrlDown=get(evd,'ControlDown');
            this.isAltDown=get(evd, 'AltDown');
            this.isMetaDown=get(evd,'MetaDown');
            this.isShiftDown=get(evd,'ShiftDown');
        end

        function nodeCollapsed(this, evd)
            uiNode  = evd.getCurrentNode;
            if this.isGroupAction && ~this.isCollapsing
                this.isCollapsing=true;
                puiNode=uiNode.getParent;
                N=puiNode.getChildCount;
                for i=0:N-1
                    this.tree.collapse(puiNode.getChildAt(i));
                end
                drawnow;
                this.isCollapsing=false;
                id=uiNode.getValue;
                this.ensureVisible(id);
            end
            
        end
        
        function setDefaultExpandedCallback(this)
            set(this.tree, 'NodeExpandedCallback', ...
                @(h,e)nodeExpanded(this, e));
        end
        
        function nodeExpanded(this, evd)
            uiNode=evd.getCurrentNode;
            if ~this.uiNodes.containsKey(uiNode.getValue)
                this.uiNodes.put(java.lang.String(uiNode.getValue), uiNode);
            end
            if ~this.tree.isLoaded(uiNode)
                children=feval(this.fncNewChildren, uiNode.getValue);
                N=length(children);
                if N==1
                    % Then we dont have an array of nodes. Create an array.
                    chnodes = children;
                    children = javaArray('com.mathworks.hg.peer.UITreeNode', 1);
                    children(1) = java(chnodes);
                end
                if N>0
                    this.tree.add(uiNode, children);
                    for i=1:N
                        this.uiNodes.put(...
                            java.lang.String(children(i).getValue), children(i));
                    end
                    this.tree.setLoaded(uiNode, true);
                end
            end
        end

        function loadChildren(this, id, force)
            if nargin<3
                force=false;
            end
            uiNode=this.uiNodes.get(id);
            if ~isempty(uiNode)
                if ~this.tree.isLoaded(uiNode) || force
                    children=feval(this.fncNewChildren, id);
                    N=length(children);
                    if N==1
                        % Then we dont have an array of nodes. Create an array.
                        chnodes = children;
                        children = javaArray('com.mathworks.hg.peer.UITreeNode', 1);
                        children(1) = java(chnodes);
                    end
                    this.tree.add(uiNode, children);
                    for i=1:N
                        this.uiNodes.put(...
                            java.lang.String(children(i).getValue), children(i));
                    end
                    this.tree.setLoaded(uiNode, true);
                end
            end
        end
        function setMultipleSelection(this,yes)
            if nargin<2
                this.tree.setMultipleSelectionEnabled(true)
            else
                this.tree.setMultipleSelectionEnabled(yes)
            end
        end
        
        function rememberNodes(this, uiNode)
            if nargin<2
                uiNode=this.root;
            end
            key=uiNode.getValue;
            this.uiNodes.put(java.lang.String(key), uiNode);
            N=uiNode.getChildCount();
            for j=1:N
                this.rememberNodes(uiNode.getChildAt(j-1));
            end
        end
        
        function ensureSelected(this, selectedIds)
            if ischar(selectedIds)
                selectedIds={selectedIds};
            end
            if ~isempty(selectedIds)
                N=length(selectedIds);
                al=java.util.ArrayList;
                for i=1:N
                    id=selectedIds{i};
                    if ~isempty(this.ensureVisible(id, false))
                        al.add(this.uiNodes.get(id));
                    end
                end
                N=al.size;
                if N>0
                    ja=javaArray('com.mathworks.hg.peer.UITreeNode', N);
                    for i=1:N
                        ja(i)=al.get(i-1);
                    end
                    this.tree.setSelectedNodes(ja);
                end
            end
        end
        
        function uiNode=ensureVisible(this, id, ...
                selectIf1MultiIf2, scroll, javaBag4Paths)
            if nargin<5
                javaBag4Paths=[];
                if nargin<4
                    scroll=true;
                    if nargin<3
                        selectIf1MultiIf2=0;
                    end
                end
            end
            if isempty(id)
                uiNode=[];
                return;
            end
            if ~this.remembered
                this.rememberNodes;
                this.remembered=true;
            end
            if ~feval(this.fncNodeExists, id)
                uiNode=[];
                warning('tree id="%s" NOT known outside of uitree...', id);
                return;
            end
            select=selectIf1MultiIf2>0;
            isMultiSelect=selectIf1MultiIf2==2;
            [inView, ~, uiNode]=SuhTree.IsInViewPort( this, id);
            if inView
                if select 
                    this.doSelection(isMultiSelect, uiNode);
                end
            else
                ok=true;
                ids=feval(this.fncGetPath, id);
                isJavaList=isa(ids, 'java.util.List');
                uiNode=this.root;
                this.expandNode(uiNode);
                if ~isJavaList
                    N=length(ids);
                    for i=1:N
                        cur=ids{i};
                        N2=uiNode.getChildCount();
                        %if N2<=0
                        %    this.loadChildren(uiNode.getValue, true);
                        %    N2=uiNode.getChildCount;
                        %end
                        ok=false;
                        for j=1:N2
                            child=uiNode.getChildAt(j-1);
                            childId=child.getValue;
                            if strcmp(childId, cur)
                                ok=true;
                                uiNode=child;
                                if i<N
                                    this.expandNode(uiNode);
                                end
                                break;
                            end
                        end
                        if ~ok
                            break;
                        end
                    end
                else
                    N=ids.size;
                    for i=0:N-1
                        cur=ids.get(i);
                        N2=uiNode.getChildCount();
                        %if N2<=0
                        %    this.loadChildren(uiNode.getValue, true);
                        %    N2=uiNode.getChildCount;
                        %end
                        ok=false;
                        for j=1:N2
                            child=uiNode.getChildAt(j-1);
                            childId=child.getValue;
                            if strcmp(childId, cur)
                                ok=true;
                                uiNode=child;
                                if i<N
                                    this.expandNode(uiNode);
                                    if uiNode.getChildCount==0
                                        this.ensureChildUiNodesExist(uiNode,true)
                                        this.expandNode(uiNode);
                                    end
                                end
                                break;
                            end
                        end
                        if ~ok
                            break;
                        end
                    end
                end
                if ~ok
                    uiNode=[];
                    warning('tree id="%s" IS known outside of uitree but NOT inside...', id);
                    return;
                elseif select 
                    this.doSelection(isMultiSelect, uiNode);
                end
                %DON'T jam it at the bottom of the window with the LEAST
                %visibility
                if nargin<4 || scroll
                    pp=uiNode.getPath;
                    tp=javax.swing.tree.TreePath(pp);
                    this.jtree.scrollPathToVisible(tp);
                    scrollUiNode=uiNode.getNextNode;
                    if ~isempty(scrollUiNode)
                        pp=scrollUiNode.getPath;
                        tp=javax.swing.tree.TreePath(pp);
                        this.jtree.scrollPathToVisible(tp);
                    end
                end
            end   
            if ~isempty(uiNode) && isjava(javaBag4Paths)
                tp=javax.swing.tree.TreePath(uiNode.getPath);
                javaBag4Paths.add(tp);
            end
            drawnow;
        end
        
        function doSelection(this, isMultiSelect, uiNode)
            if isMultiSelect
                tp=javax.swing.tree.TreePath(uiNode.getPath);
                this.jtree.addSelectionPath(tp);                
            else
                this.tree.setSelectedNode(uiNode);
            end
        end
        
        function expandNode(this, uiNode)
            try
                this.ensureChildUiNodesExist(uiNode);
                if ~uiNode.isLeafNode()
                    this.tree.expand(uiNode);
                    drawnow;
                end
            catch ex
                ex.getReport
            end
        end
        
        function [ok, uiNode]=ensureChildUiNodesExist(this, uiNode, force)
            if nargin<3
                force=false;
            end
            if ischar(uiNode)
                uiNode=this.uiNodes.get(java.lang.String(uiNode));
                if isempty(uiNode)
                    ok=false;
                    return;
                end
            end
            ok=true;
            uiN=uiNode.getChildCount();
            id=uiNode.getValue;
            idChildren=feval(this.fncGetChildren, id);
            nonUiN=length(idChildren);
            if uiN~=nonUiN || force
                try
                    [uiChildren, sortedIdChildren]=feval(this.fncNewChildren, id);
                catch 
                    uiChildren=feval(this.fncNewChildren, id);
                    sortedIdChildren=idChildren;
                end
                nUiChildren=length(uiChildren);
                if uiNode.isLeafNode
                    uiNode.setLeafNode(false);
                end
                %if uiN==0
                %    disp('herp');
                %end
                if uiN==0 || nonUiN<uiN || force
                    uiNode.removeAllChildren;
                    for i=1:nUiChildren
                        uiNode.add(uiChildren(i));
                        this.uiNodes.put(java.lang.String( ...
                            sortedIdChildren{i}), uiChildren(i));
                    end
                else
                    orderChanged=false;
                    if nUiChildren>0
                        for i=1:nonUiN
                            if i>nUiChildren
                                break;
                            end
                            idChild=java.lang.String(sortedIdChildren{i});
                            %if strcmp(idChild, 'gate:ID2123123387')
                            %    disp('herp')
                            %end
                            if ~this.uiNodes.containsKey(idChild)
                                uiNode.add(uiChildren(i));
                                this.uiNodes.put(idChild, uiChildren(i));
                            elseif ~orderChanged
                                uiChild2=uiNode.getChildAt(i-1);
                                orderChanged=~isequal( ...
                                    uiChild2.getValue, idChild);
                            end
                        end
                    end
                    if orderChanged
                        newTotal=uiNode.getChildCount;
                        if newTotal==nUiChildren
                            for i=1:newTotal
                                uiChild=uiChildren(i);
                                uiChild2=uiNode.getChildAt(i-1);
                                if ~isequal(uiChild2.getValue, ...
                                        uiChild.getValue)
                                    uiChild2.setName(uiChild.getName);
                                    uiChild2.setValue(uiChild.getValue);
                                    if uiChild2.getChildCount>0
                                        fprintf(['"%s" overlays ' ...
                                            '"%s" position\n... so ' ...
                                            'remove %d children\n'], ...
                                            Html.Remove(uiChild.getName), ...
                                            Html.Remove(uiChild2.getName), ...
                                            uiChild2.getChildCount)
                                        uiChild2.removeAllChildren;
                                        this.tree.setLoaded(uiChild2, false);
                                    end
                                end
                                this.uiNodes.put(java.lang.String(...
                                    uiChild.getValue), uiChild2);
                            end
                        end
                    end
                end
                this.tree.reloadNode(uiNode);
                drawnow;
                this.tree.setLoaded(uiNode,true);
            end
        end
        
        function refreshNode(this, uiNode, htmlOrText)
            uiNode.setName(htmlOrText);
            %was=this.IsExpanded(this.jtree, uiNode);
            this.tree.reloadNode(uiNode);
        end
        
        function stylize(this, fontName)
            if nargin<2
                fontName='Arial';
            end
            jt=this.jtree;
            jt.setRowHeight(0);
            f=jt.getFont;
            jt.setFont( java.awt.Font(fontName, f.getStyle, f.getSize-1));
            jt.setBorder(javax.swing.BorderFactory.createEmptyBorder (12,12,12,12));
            jt.setRowHeight(0);
        end
        
    end
    
    methods(Static)
        function this=New(root, fncNodeSelected, fncGetPath, fncNodeExists, ...
                fncNewChildren, fncGetChildren, visible)
            if nargin<7
                visible=true;
            end
            [uit, container] = uitree('v0','Root',  root, ...
                'SelectionChangeFcn',fncNodeSelected);
            uit.Visible=visible;
            this=SuhTree(uit, container, root, fncGetPath, ...
                fncNodeExists, fncNewChildren, fncGetChildren);
        end
        
        function [ok, selected, uiNode, row]=IsInViewPort(this, id)
            ok=false;
            row=0;
            selected=false;
            uiNode=this.uiNodes.get(java.lang.String(id));
            if isempty(uiNode)
                return;
            end
            tr=this.jtree;
            vp=javaObjectEDT(tr.getParent);
            vr=javaObjectEDT(vp.getViewRect);
            firstRow=tr.getClosestRowForLocation( vr.x, vr.y);
        	lastRow=tr.getClosestRowForLocation(vr.x, vr.y + vr.height);   
            for i=firstRow+1:lastRow-1
                path=javaObjectEDT(tr.getPathForRow(i));
                uiNode=path.getLastPathComponent;
                id1=uiNode.getValue;
                if strcmp(id, id1)
                    ok=true;
                    selected=tr.isRowSelected(i);
                    row=i;
                    return;
                end
            end
        end
        
        function uiNode=GetClickedNode(eventData)
            treePath=eventData.getSource.getPathForLocation(...
                eventData.getX, eventData.getY);
            if ~isempty(treePath)
                uiNode=treePath.getLastPathComponent;
            else
                uiNode=[];
            end
        end
        
        function uiNode=ClickedBottomRight(eventData)
            uiNode=[];
            clickX = eventData.getX;
            clickY = eventData.getY;
            jtree_ = eventData.getSource;
            treePath = jtree_.getPathForLocation(clickX, clickY);
            if ~isempty(treePath)
                testPath = jtree_.getPathForLocation(clickX+62, clickY);
                if isempty(testPath)
                    testPath=jtree_.getPathForLocation(clickX, clickY+40);
                    if isempty(testPath) || testPath ~= treePath
                        uiNode = treePath.getLastPathComponent;
                        %disp(['clicked BOTTOM RIGHT of ' char(uiNode.getName)]);
                    end
                end
            end
        end
        
        function yes=IsExpanded(jtree, uiNode)
            path=uiNode.getPath;
            tp=javax.swing.tree.TreePath(path);
            yes=jtree.isExpanded(tp);
        end

        function expanded=GetExpanded(jtree, uiNode, expanded)
            if nargin<3
                expanded=java.util.ArrayList;
            end
            if SuhTree.IsExpanded(jtree, uiNode)
                expanded.add(uiNode);
            end
            N=uiNode.getChildCount();
            for i=1:N
                uiChild=uiNode.getChildAt(i-1);
                expanded=SuhTree.GetExpanded(jtree, uiChild, expanded);
            end
        end
        
        function unexpanded=GetUnexpanded(jtree, root, uiNode, unexpanded)
            if nargin<4
                unexpanded=java.util.ArrayList;
            end
            N=uiNode.getChildCount();
            if ~uiNode.equals(root) && ...
                    (N==0 || ~SuhTree.IsExpanded(jtree, uiNode))
                unexpanded.add(uiNode);
            else
                for i=1:N
                    uiChild=uiNode.getChildAt(i-1);
                    unexpanded=SuhTree.GetUnexpanded(jtree, root, uiChild, unexpanded);
                end
            end
        end
    end
end