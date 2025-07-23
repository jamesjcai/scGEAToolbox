function [answer, cancelled]=inputDlg(msg,title,varargin)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

if nargin==1
    title='Input required';
end
[msgType, jsa,default,~,~]=getMsgTypeAndOptions(...
    javax.swing.JOptionPane.QUESTION_MESSAGE, varargin);
[msg, where, property, properties, ~, myIcon, javaWin, tip,~,modal]...
    =decodeMsg(msg, default);
hasProp=~isempty(properties) && ~isempty(property);
inputValue='';
if hasProp
    inputValue=properties.get(property);
end
if isempty(inputValue) && ~isempty(jsa)
    inputValue=char(jsa(1));
end
if isempty(myIcon)
    if msgType==0
        myIcon='error.png';
    elseif msgType==1
        myIcon = 'facs.gif';
    elseif msgType==2
        myIcon='warning.png';
    else
        myIcon='question.png';
    end
end
pane=javaObjectEDT('javax.swing.JOptionPane', msg, msgType);
pane.setWantsInput(true);
pane.setInitialSelectionValue(inputValue);
pane.selectInitialValue();
pane.setIcon(Gui.Icon(myIcon));
pane.setOptionType(javax.swing.JOptionPane.OK_CANCEL_OPTION);
if ~isempty(tip)
    pane.setToolTipText(tip);
end
MatBasics.RunLater(@(h,e)onTop(),.5);
PopUp.Pane(pane, title,where, javaWin, modal);
answer=pane.getInputValue;
cancelled=strcmp(answer,'uninitializedValue');
if cancelled
    answer='';
else
    if hasProp
        properties.set(property, answer);
    end
end
    function onTop
        setAlwaysOnTopTimer( Gui.WindowAncestor(pane), 5, true);
    end
end
