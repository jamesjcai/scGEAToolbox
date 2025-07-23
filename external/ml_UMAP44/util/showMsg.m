function pu=showMsg(msg, title, where, modal, showBusy, pauseSecs, addClose)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

if nargin<7
    addClose=true;
    if nargin<6
        pauseSecs=0;
        if nargin<5
            showBusy=true;
            if nargin<4
                modal=true;
                if nargin<3
                    where='center';
                    if nargin<2
                        title='Note....';
                    end
                end
            end
        end
    end
end
pu=PopUp(msg, where, title, showBusy, [],[],modal);
if addClose
    pu.addCloseBtn;
end
if modal
    pu.dlg.setAlwaysOnTop(true);
    Gui.SetJavaVisible(pu.dlg);
    
elseif pauseSecs>0
    PopUp.TimedClose(pu.dlg, pauseSecs);
end
