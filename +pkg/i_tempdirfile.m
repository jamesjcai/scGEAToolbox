function [a, b] = i_tempdirfile(subdirname, ext)
if nargin < 2
    ext = 'txt';
end
if nargin < 1
    subdirname = [];
end

if isempty(subdirname)
    if ~isMATLABReleaseOlderThan('R2025a')
        subdirname = ['pid_' num2str(matlabProcessID())];
    else
        subdirname = ['pid_' num2str(feature('getpid'))];
    end
end
a = fullfile(tempdir, subdirname);
if ~exist(a, "dir")
    mkdir(a);
end

if nargout > 1
    tstr = matlab.lang.makeValidName(string(datetime));
    b = fullfile(a, strcat(tstr, ".", string(ext)));
end

end
