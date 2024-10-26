function answer=ask(question, defaultAnswer, title, where, cancelToo)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab 
%   License: BSD 3 clause

if nargin<5
    cancelToo=true;
    if nargin<4
        where='center';
        if nargin<3
            title=[];
            if nargin<2
                defaultAnswer=1;
            end
        end
    end
end
if isempty(title)
    title='Please confirm....';
end
pu=PopUp(question, where, title, false, [],[], true);
pu.addYesNo(cancelToo, defaultAnswer);
pu.dlg.setAlwaysOnTop(true);
pu.dlg.setModal(true);
Gui.SetJavaVisible(pu.dlg);
answer=pu.answer;
