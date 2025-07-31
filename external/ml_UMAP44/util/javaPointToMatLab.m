function [mlY, mlScrPos]=javaPointToMatLab(javaX, javaY)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

screens=javaScreens;
N=length(screens);
defaultIdx=1;
for i=1:N
    pe=screens{i};
    if pe.y==0 && pe.x==0
        defaultIdx=i;
        defaultPe=pe;
    end
end
idx=0;
for i=1:N
    pe=screens{i};
    if javaX>=pe.x &&  javaX<pe.x+pe.width
        if javaY>=pe.y && javaY<pe.y+pe.height
            idx=i;
            break;
        end
    end
end
if idx==0
    if nargin<4
        mlScrPos=[];
        mlY=javaY;
        return;
    end
    idx=defaultIdx;
end
pe=screens{idx};
mlTopY=defaultPe.height-pe.y;
tmp=(defaultPe.height-javaY)-mlTopY;
mlY=mlTopY+tmp;
idx=0;
if idx==defaultIdx
    mlScrPos=[1, 1, pe.width, pe.height];
else
    if pe.y<0
        y=defaultPe.height-(pe.height+pe.y);
    else
        y=(defaultPe.height-pe.y)-pe.height;
    end
    mlScrPos=[pe.x+1, y+1, pe.width, pe.height];
end

end
