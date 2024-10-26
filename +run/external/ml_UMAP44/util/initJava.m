function ok=initJava
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

ok=false;
try
    edu.stanford.facs.swing.StochasticGradientDescent.EPOCH_REPORTS;
catch
    jar=fullfile(fileparts(mfilename('fullpath')), 'suh.jar');
    javaaddpath(jar);
end
try
    edu.stanford.facs.swing.StochasticGradientDescent.EPOCH_REPORTS;
    ok=true;
catch
end