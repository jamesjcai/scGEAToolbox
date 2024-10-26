function val = clip(val)
%CLIP Standard clamping of a value into a fixed range (in this case -4 to
% 4). For speed reasons, this function is no longer used by the MATLAB
% algorithm.
%
% val = CLIP(val)
% 
% Parameters
% ----------
% val: double
%     The value to be clamped.
% 
% Returns
% -------
% The clamped value, now fixed to be in the range -4 to 4.

%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

 val = min(4, val);
 val = max(-4, val);