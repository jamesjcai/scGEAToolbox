function [a, b] = i_tempdirfile(subdirname, ext)
if nargin < 2
    ext = 'txt';
end
if nargin < 1
    subdirname = [];
end

if ~isempty(subdirname)
    a = fullfile(tempdir, subdirname);
    if ~exist(a, "dir")
        mkdir(a);
    end
else
    a = tempdir;
end

if nargout > 1
    tstr = matlab.lang.makeValidName(string(datetime));
    b = fullfile(a, strcat(tstr, ".", string(ext)));
end

end