function [jd,pane]=msgWarning(txt, pause, where, title)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

if nargin<4
    title='Warning...';
    if nargin<3
        where='center';
        if nargin<2
            pause=9;
        end
    end
end
[jd, pane]=msg(txt, pause, where, title, 'warning.png');
try
    if isstruct(txt)
        warning(txt.msg);
    else
        warning(txt);
    end
catch ex
    ex.getReport
end
Gui.Boing();
end
