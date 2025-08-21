function addDLLLocationToSystemPathIfNecessary()
%ADDDLLLOCATIONTOSYSTEMPATHIFNECESSARY
%  On Windows, adds the location of the necessary DLLs distributed with this
%  support package to the system path.  Without this, the mex functions
%  distributed with this support package do not function.
%
%  For Linux/Mac platforms, this is a no-op.
%

% Copyright 2023 The MathWorks, Inc.

if ~ispc()
    return;
end

persistent initd
mlock
if isempty(initd)
    pathEntryToAdd = fullfile(...
        matlabshared.supportpkg.getSupportPackageRoot, "bin", computer("arch"));
    existingPath = getenv("PATH");
    existingPathEntries = strsplit(existingPath, pathsep);
    if ~any(strcmpi(existingPathEntries, pathEntryToAdd))
        newPath = strcat(existingPath, pathsep, pathEntryToAdd);
        setenv("PATH", newPath);
    end
    initd = true;
end

