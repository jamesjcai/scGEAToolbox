function p = i_normalizepath(p)
%NORMALIZEPATH Convert a path containing ".." or "." into a canonical path
%
%   p = normalizePath(p)
%
%   Example:
%   p = normalizePath("D:\GitHub\scGEAToolbox_dev\+run\..\external\py_cellbender\script.py")

    if isstring(p) || ischar(p)
        p = char(p);
    else
        error('Input must be a string or char array.');
    end

    try
        % Use Java canonical path resolution
        f = java.io.File(p);
        p = char(f.getCanonicalPath());
    catch
        % Fallback (if Java fails for some reason)
        p = char(java.io.File(p).getAbsolutePath());
    end
end
