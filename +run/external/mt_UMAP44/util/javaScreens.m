function screens=javaScreens
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    screens={};
    ge=java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment;
    physicalScreens=ge.getScreenDevices;
    N=length(physicalScreens);
    for i=1:N
        pe=physicalScreens(i).getDefaultConfiguration.getBounds;
        screens{end+1}=pe;
    end
end
