function mustBeThreshold(thr)
if ~(matlab.internal.math.partialMatch(thr,"none") || (isscalar(thr) && isnumeric(thr) && isreal(thr) && ~isnan(thr)))
    error(message('quantum:QuantumState:invalidThreshold'))
end
end