function [yes, cancelled, raIfRemembered, jd]...
    =askYesOrNo(theMsg, title, where, ...
    defaultIsYes, rememberId, property)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab 
%   License: BSD 3 clause

if nargin<4
    defaultIsYes=true;
    if nargin<3
        where=[];
        if nargin<2
            title=[];
        end
    end
end
if isempty(title)
    title= 'Please confirm...';
end
if nargin>2
    if isempty(where)
        where='center';
    end
    if ~isstruct(theMsg)
        m.msg=theMsg;
        m.where=where;
        if nargin>4 
            if ~isempty(rememberId)
                m.remember=rememberId;
            end
            if nargin>5
                m.property=property;
            end
        end
        theMsg=m;
    else
        theMsg.where=where;
        if nargin>4 && ~isempty(rememberId)
            theMsg.remember=rememberId;
        end
        if nargin>5
            theMsg.property=property;
        end
    end
end
if defaultIsYes
    dflt='Yes';
else
    dflt='No';
end
if nargout>1
    [~,yes,cancelled, raIfRemembered, jd]=questDlg(theMsg, ...
        title, 'Yes', 'No', 'Cancel', dflt);
else
    [~,yes,cancelled, raIfRemembered, jd]=questDlg(theMsg, ...
        title, 'Yes', 'No', dflt);
end
end