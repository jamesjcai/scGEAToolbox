function mustBeRotationThreshold(thr)
%

% Copyright 2023-2024 The MathWorks, Inc.

if matlab.internal.math.partialMatch(thr,"none") || matlab.internal.math.partialMatch(thr,"auto")
    return
elseif (isscalar(thr) && isnumeric(thr) && isreal(thr) && thr>=0)
    return
else
    error(message('quantum:CompositeGate:invalidRotationThreshold'))
end
end
