function [x, y]=clockwise(xy)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

x=xy(:,1);
y=xy(:,2);
cx = mean(x);
cy = mean(y);
a = atan2(y - cy, x - cx);
[~, order] = sort(a);
x = x(order);
y = y(order);
%assert(all(ismember([x y],xy,'rows')));
end
