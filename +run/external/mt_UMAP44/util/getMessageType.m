function [type, found]=getMessageType(arg, dflt)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

if nargin<2
        dflt=javax.swing.JOptionPane.PLAIN_MESSAGE;
end
found=true;
if strcmpi(arg, 'error')
        type=javax.swing.JOptionPane.ERROR_MESSAGE;
    elseif strcmpi(arg, 'question')
        type=javax.swing.JOptionPane.QUESTION_MESSAGE;
    elseif strcmpi(arg, 'warning') || strcmpi(arg, 'warn')
        type=javax.swing.JOptionPane.WARNING_MESSAGE;
    elseif strcmpi(arg, 'information')
        type=javax.swing.JOptionPane.INFORMATION_MESSAGE;
    elseif strcmpi(arg, 'plain')
        type=javax.swing.JOptionPane.PLAIN_MESSAGE;
else
        found=false;
        type=dflt;
end
end