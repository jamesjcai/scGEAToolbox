function filepath = defaultFilepath(filename)
% Internal use only.

% Copyright 2023 The MathWorks, Inc.
if ispc
    homepath = getenv("USERPROFILE");
elseif ismac
    homepath = getenv("HOME");
else
    assert(isunix)
    homepath = getenv("HOME");
end
filepath = fullfile(homepath, ".qiskit", filename);
end
