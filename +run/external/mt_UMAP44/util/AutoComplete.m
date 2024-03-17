function [jSearch, researchFnc, jCombo, hSearch, hCombo]=...
    AutoComplete(findFncOrItemsToFind, panel, position, label, ...
     doneButton, escButton, chbContains, minimumEntries )
%Generic auto complete widget derived from Yair Altman's example found at 
% http://undocumentedmatlab.com/blog/auto-completion-widget.
% Enhancements here make the widget more application-independent and 
% more OS independent (e.g. Compatible with Mac's aqua L&F).  
% Yair's example seems specific to asset lookup, to MS Windows

%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

if nargin<8
    minimumEntries=0;
    if nargin<7
        chbContains=[];
        if nargin<6
            escButton=[];
            if nargin<5
                doneButton=[];
                if nargin<4
                    label=[];
                    if nargin<3
                        position=[];
                        if nargin<2
                            panel=uipanel;%('position', [25, 25, 100, 20]);
                        end
                    end
                    
                end
            end
        end
    end
end

if iscell(findFncOrItemsToFind)
    items=findFncOrItemsToFind;
    findFncOrItemsToFind=@findItems;
end
    function l=findItems(searchText)
        xx=regexpi(items, searchText);
        l=items(~cellfun(@isempty, xx)); 
    end
if isempty(position)
    position=[1, 1, 0, 25];
end
if isempty(label)
    if ispc
        label='Enter search value(s) ...';
    else
        label='Enter search value...';
    end
end
isJavaContainer=isa(panel, 'javax.swing.JPanel');
jCombo=javaObjectEDT('javax.swing.JComboBox');
Y=position(1);
X=position(2);
if position(3)<1
    position(3)=length(label)*7;
end
W=position(3);
H=position(4);
fgNotFound=java.awt.Color(.6, .2, 1);
fgNotFoundStyle=java.awt.Font.BOLD+java.awt.Font.ITALIC;    
if ispc
    bgNotFound=java.awt.Color(243/255, 252/255, .8);
    popupMargin=1;
else
    bgNotFound=java.awt.Color(243/255, 252/255, .7);
    popupMargin=20;
end
mi=javaObjectEDT(javax.swing.plaf.metal.MetalComboBoxUI);
mi.installUI(jCombo);
if ~isJavaContainer
      [jCombo, hCombo]=javacomponent(jCombo, ...
        [X+popupMargin Y W-popupMargin H],panel); %#ok<*JAVCM> 
end
jCombo.getComponent(2).setVisible(false);
jSearch=javaObjectEDT(com.mathworks.widgets.SearchTextField(label)); %#ok<*JAPIMATHWORKS> 
jSearchComponent=jSearch.getComponent;
jSearchButton=handle(jSearchComponent.getComponent(1), 'CallbackProperties');
set(jSearchButton, 'ActionPerformedCallback',  {@searchKeyCb, jCombo, jSearch});
%jSearchButton.setEnabled(true);
hCombo2=handle(jCombo, 'CallbackProperties');

if ispc
    hSearchField = handle(jSearchComponent.getComponent(0), ...
        'CallbackProperties');
    set(hCombo2, 'KeyPressedCallback', {@comboKeyCb, jCombo, jSearch});
    set(hCombo2, 'ActionPerformedCallback', @(src,evd)selectComboCb(src, ...
        evd, jCombo,jSearch));    
    comboOffsetY=4;
    comboOffsetWidth=26;
else
    hSearchField = handle(jSearchComponent, 'CallbackProperties');
    set(jCombo, 'ActionPerformedCallback', @(src,evd)selectComboCb(src, ...
        evd, jCombo,jSearch));    
    comboOffsetY=2;
    comboOffsetWidth=20;
end
set(hSearchField, 'KeyPressedCallback', {@searchKeyCb, jCombo, jSearch});
jSearchComponent.requestFocus
lastSearchText = '';  
if isJavaContainer
    hSearch=[];
    hCombo=[];
    jp2=javaObjectEDT('javax.swing.JPanel');
    jp2.add(jSearchComponent);
    jSearchComponent.setPreferredSize(...
        java.awt.Dimension(W, H));
    if ispc
        layout=javaObjectEDT(...
            edu.stanford.facs.swing.XYLayout.layoutAutoComplete(panel,...
            jp2, jCombo, X, Y, W, H, 0, Y+comboOffsetY));
    else
        layout=javaObjectEDT(...
            edu.stanford.facs.swing.XYLayout.layoutAutoComplete(panel,...
            jp2, jCombo, X, Y, W, H, comboOffsetWidth, Y+comboOffsetY));
    end
else
    [~, hSearch] = javacomponent(jSearchComponent,[X,Y,W,H],panel);
end
researchFnc=@(resetLastSearched, kc, textOverride)searchCb(...
    resetLastSearched, false, kc, jCombo, jSearch, textOverride);

    function comboKeyCb(~, eventData,jCombo, jSearch)
        lastSearchText=jCombo.getSelectedItem;
        kc=eventData.getKeyCode;
        if kc==10 %ENTER pressed
            lastSearchText=jCombo.getSelectedItem;
            jSearch.setSearchText(lastSearchText);
            jCombo.hidePopup;
            if ispc
                jCombo.setVisible(false);
            end
            jSearch.requestFocus;
        elseif kc==27
            jCombo.hidePopup;
            if ispc
                jCombo.setVisible(false);
            end
        end
    end

    function selectComboCb(~,~,jCombo, jSearch)
        lastSearchText=jCombo.getSelectedItem;
        jSearch.setSearchText(lastSearchText)
        hSearchField.requestFocus;
        if ispc
            jCombo.setVisible(false);
        end
        
    end

    function searchKeyCb(~, eventData, jCombo, jSearch)
        try
            kc=eventData.getKeyCode;
            searchButtonClicked=false;
        catch 
            searchButtonClicked=true;
            kc=0;
        end
        searchCb(false, searchButtonClicked, kc, jCombo, jSearch);
    end

    function searchCb(resetLastSearch, searchButtonClicked, ...
            kc, jCombo, jSearch, textOverride)
        pressedDown=kc==40;
        pressedEsc=kc==27;
        pressedEnter=kc==10;
        if isempty(lastSearchText) || resetLastSearch
            lastSearchText = '';
        end
        drawnow;
        if nargin<6|| isempty(textOverride)
            enteredText=char(jSearch.getSearchText);
            N=length(enteredText);
            if minimumEntries~=0  && N<=minimumEntries && ~pressedDown
                return;
            end
        else
            enteredText=textOverride;
            N=length(enteredText);
        end
        empty=isempty(N==0);        
        if ~isempty(chbContains) && ((islogical(chbContains) ...
                && chbContains)  || (~islogical(chbContains) && ...
                get(chbContains, 'value')==1))
            searchText=['.*' strrep(enteredText, '*', '.*') '.*'];  % turn into a valid regexp
        else
            searchText=['^' strrep(enteredText, '*', '.*')];  % turn into a valid regexp
        end
        
        searchText=strrep(searchText, '(', '\(');
        searchText=strrep(searchText, ')', '\)');
        searchText=strrep(searchText, '+', '\+');
        disp(searchText);
        searchText=regexprep(searchText, '<[^>]>', '');
        if strcmpi(searchText, lastSearchText) && ~isempty(searchText)
            done=true;
            if pressedDown
                if strcmp(searchText, '^') && ~jCombo.isPopupVisible
                    searchText='^.*';
                    done=false;
                end                
            end            
            if jCombo.getItemCount>0
                if pressedEnter
                    if jCombo.isPopupVisible
                        lastSearchText=jCombo.getSelectedItem;
                        jSearch.setSearchText(lastSearchText);
                    elseif ~isempty(doneButton)
                        doneButton.doClick;
                        return;
                    end
                    jCombo.hidePopup;
                    if ispc
                        jCombo.setVisible(false);
                    end
                else
                    if ispc
                        jCombo.setVisible(true);
                        drawnow;
                    end
                    jCombo.showPopup;
                    if pressedDown
                        jCombo.requestFocus
                    end
                end
            end
            if done
                if pressedEnter
                    if ~isempty(doneButton)
                        doneButton.doClick;
                    end
                elseif pressedEsc
                    if jCombo.isPopupVisible
                        jCombo.hidePopup;
                        if ispc
                            jCombo.setVisible(false);
                        end
                        return;
                    elseif ~isempty(escButton)
                        escButton.doClick;
                        return;
                    end
                end
                return;  % maybe just clicked an arrow key or Home/End - no need to refresh the popup panel
            end
        elseif pressedEsc
            if jCombo.isPopupVisible
                jCombo.hidePopup;
                if ispc
                    jCombo.setVisible(false);
                end
                return;
            elseif ~isempty(escButton)
                escButton.doClick;
                return;
            end
        elseif pressedEnter
            if ~isempty(doneButton)
                doneButton.doClick;
                return;
            end
        elseif (pressedDown || searchButtonClicked) && strcmp(searchText, '^') && ~jCombo.isPopupVisible
            searchText='^.*';
        end
        if nargin<6
            lastSearchText = searchText;
        end
        vector=findFncOrItemsToFind(searchText);
        if iscell(vector)
            N=length(vector);
        else
            N=vector.size;
        end
        if N>0
            if ~isa(vector, 'java.util.Vector')
                v=java.util.Vector;
                v.addAll(vector);
                vector=v;
            end
            jCombo.setModel(javax.swing.DefaultComboBoxModel(vector));
            if kc<0
                return;
            end
            if ispc
                jCombo.setVisible(true);
                drawnow;
            end
            if ismac
                if isJavaContainer
                    if N>jCombo.getMaximumRowCount
                        %disp(layout.constraints);
                        layout.setConstraint(jCombo, 20, 5);
                        %disp('NOW');
                        %disp(layout.constraints);
                    else
                        layout.setConstraint(jCombo, 20, 24);
                    end
                    drawnow;
                else
                    pos=get(hCombo, 'position');
                    pos2=get(hSearch, 'position');
                    if N>jCombo.getMaximumRowCount
                        set(hCombo, 'position', [pos(1) pos2(2) pos(3) pos(4)]);
                    else
                        set(hCombo, 'position', [pos(1) pos2(2)-19 pos(3) pos(4)]);
                        drawnow;
                    end
                end
            end
            stylize(java.awt.Font.PLAIN, java.awt.Color.white, ...
                    java.awt.Color.blue);
            jCombo.showPopup;    
            if pressedDown
                jCombo.requestFocus
            end
        else
            jCombo.setModel(javax.swing.DefaultComboBoxModel(java.util.Vector));
            if empty
                stylize(java.awt.Font.PLAIN, java.awt.Color.white, ...
                    java.awt.Color.blue);
            else
                stylize(fgNotFoundStyle, bgNotFound, fgNotFound);
            end
            jCombo.hidePopup;
            if ispc
                jCombo.setVisible(false);
            end
        end
    end  

    function stylize(style, bg, fg)        
        hSearchField.setFont(hSearchField.getFont.deriveFont(uint8(style)));
        hSearchField.setBackground(bg);
        hSearchField.setForeground(fg);
    end
end