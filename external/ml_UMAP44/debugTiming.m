function debugTiming(description)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

DEBUGGING=false;
if DEBUGGING && ~isdeployed
    msg([description ' ' num2str(toc)], 40, ...
        'south west+', 'Timing test', 'none')
end
end