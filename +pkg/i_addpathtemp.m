function cleanupobj = i_addpathtemp(varargin)
% I_ADDPATHTEMP - Add folders to the search path for the caller's lifetime
%
%   C = PKG.I_ADDPATHTEMP(FOLDER1, FOLDER2, ...) adds each named folder that
%   is not already on the MATLAB search path and returns an ONCLEANUP object.
%   When C goes out of scope -- normally when the calling function returns,
%   errors, or is interrupted -- the folders that were added are removed
%   again, so the caller leaves the path as it found it. Folders that were
%   already on the path are never removed.
%
%   Use this for third-party folders whose file names are generic enough to
%   shadow toolbox or user code (for example external/ml_umap45/util, which
%   defines Args.m, File.m, Map.m, Plots.m and String.m).
%
%   Assign the output to a variable that lives as long as the folders are
%   needed. Calling PKG.I_ADDPATHTEMP(...) without capturing the output
%   restores the path immediately and so has no useful effect.
%
%   Example:
%   c = pkg.i_addpathtemp(pth, fullfile(pth, 'util'));   %#ok<NASGU>
%
%   See also ADDPATH, ONCLEANUP.

folders = cellfun(@pkg.i_normalizepath, varargin, UniformOutput=false);

% Only take responsibility for folders this call actually adds; one that was
% already on the path belongs to whoever put it there.
onpath = strsplit(path, pathsep);
if ispc
    isnew = ~ismember(lower(folders), lower(onpath));
else
    isnew = ~ismember(folders, onpath);
end
folders = folders(isnew);

if ~isempty(folders)
    addpath(folders{:});
end
cleanupobj = onCleanup(@() in_restorepath(folders));
end

function in_restorepath(folders)
if isempty(folders)
    return;
end
% The search path can legitimately have changed while the caller ran, so a
% folder that is no longer there is not worth a warning.
w = warning('off', 'MATLAB:rmpath:DirNotFound');
restorewarning = onCleanup(@() warning(w));
try
    rmpath(folders{:});
catch
    % Failing to shrink the path is not worth surfacing: the caller's result
    % is already computed, and the only cost is a stale path entry.
end
end
