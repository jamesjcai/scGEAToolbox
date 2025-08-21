function str = versionString()
%

% Copyright 2023 The MathWorks, Inc.

s = matlabshared.supportpkg.getInstalled;
if isempty(s)
    % No support packages installed - send MATLAB version instead
    s = ver('MATLAB'); %#ok<VERMATLAB>
    spkgVersion = s.Version;
else
    names = string({s.Name});
    versions = {s.InstalledVersion};
    spkgVersion = versions{"MATLAB Support Package for Quantum Computing" == names};
end
str = "MATLAB_SPKG_QUANTUM/" + spkgVersion;
