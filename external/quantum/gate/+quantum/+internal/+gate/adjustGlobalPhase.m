function [phase, err] = adjustGlobalPhase(actual, expected)
%ADJUSTGLOBALPHASE Adjust global phase to best match expected value
%
%   FOR INTERNAL USE ONLY -- This feature is intentionally undocumented.
%   Its behavior may change, or it may be removed in a future release.
%
%   phase = adjustGlobalPhase(actual, expected) returns a global phase to
%   be applied to actual to make it match expected as closely as possible.
%   Inputs actual and expected can be vectors or matrices of complex
%   numbers.
%
%   Output phase is a scalar with absolute value 1 that minimizes
%              norm(phase*actual - expected, 'fro')
%
%   [phase, err] = adjustGlobalPhase(actual, expected) additionally returns
%   the error err = norm(phase*actual - expected, 'fro').
%
%   In quantum computing, the global phase is not measurable and therefore
%   we often want to compare two states without taking interest in the
%   global phase.
%
%   See also quantum.gate.QuantumState

%   Copyright 2023-2025 The MathWorks, Inc.

assert(isequal(size(actual), size(expected)))

phase = sign(actual(:)' * expected(:));

if phase == 0
    phase = 1;
end

err = norm(expected - phase*actual, 'fro');
end