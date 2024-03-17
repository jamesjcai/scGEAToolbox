function bins=probability_means_weights_ptrs(data)
% DEPRECATED .... use util/SuhProbabilityBins(data)

%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab 
%   License: BSD 3 clause

MIN_BINS=8192;
MIN_EVENTS_PER_BIN=4;
MAX_EVENTS_PER_BIN=34;
N=size(data, 1);
events_per_bin=floor(2*log(N));
number_of_bins=floor(N/events_per_bin);
if number_of_bins<MIN_BINS
    events_per_bin=floor(N/MIN_BINS);
end
if events_per_bin<MIN_EVENTS_PER_BIN
    events_per_bin=MIN_EVENTS_PER_BIN;
elseif events_per_bin>MAX_EVENTS_PER_BIN
    events_per_bin=MAX_EVENTS_PER_BIN;
end
if number_of_bins>2^14 %16384
    events_per_bin=MAX_EVENTS_PER_BIN;
end
[bins.means, bins.ptrs, ~, bins.weights]=...
    probability_bin(data, data, events_per_bin, false);

end