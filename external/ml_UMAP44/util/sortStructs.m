function [out, N]=sortStructs(in, field, direction)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

out=in;
N=length(in);
if N>0
    value=out(1).(field);
    if isnumeric(value)
        S=[out(:).(field)];
        [~,si]=sort(S, direction);
    elseif ischar(value)
        S={out(:).(field)};
        [~,si]=sort(S);
        if isequal(direction, 'descend')
            si=flip(si);
        end
    end
    out=out(si);
end
end
