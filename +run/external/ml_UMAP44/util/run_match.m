function [match, matchTable, trainingQfTreeMatch, ...
    trainingQfTree, testQfTreeMatch,  testQfTree]...
    =run_match(varargin)
%   Wrapper for SuhMatch.Run so that the command line syntax can be used 
%   without extra characters such as ( or ) or , or '

%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab 
%   License: BSD 3 clause
    [match, matchTable, trainingQfTreeMatch, ...
        trainingQfTree, testQfTreeMatch, testQfTree]...
        =SuhMatch.Run(varargin{:});
end