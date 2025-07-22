function [d, resized]=resizeJavaPrefs(app, J, H)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    d=J.getPreferredSize;
    if app.toolBarFactor>0
        d=java.awt.Dimension(d.width*.45, d.height*.6);
        J.setPreferredSize(d);
        resized=true;
    else
        resized=false;
    end
end
