function [jd,pane]=msgError(txt, pause, where, title)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab 
%   License: BSD 3 clause

if nargin<4
    title='Error...';
    if nargin<3
        where='center';
        if nargin<2
            pause=0;
        end
    end
end
if strcmpi('modal', pause)
    Gui.Splat;
    [jd, pane]=msg(struct('modal',  true, 'msg', txt),...
        0, where, title, 'error.png');
else
    if ~isnumeric(pause)
        pause=8;
    end
    [jd, pane]=msg(txt, pause, where, title, 'error.png');
    Gui.Splat(jd);
end
end
