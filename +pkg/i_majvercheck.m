function [needupdate, v1local, v2web, im] = i_majvercheck(needed)

if nargin<1, needed = [true true true true]; end
% major version update check
needupdate = false;
v1local = [];
v2web = [];
im = [];
try
    if needed(2)
        v1local = pkg.i_get_versionnum;
    end
    if needed(3)

        instURL = 'https://api.github.com/repos/jamesjcai/scGEAToolbox/releases/latest';
        instRes = webread(instURL);
        v2web = instRes.tag_name(2:end);

    end
    % Only flag an update when the released version is strictly newer than
    % the local one. A local version ahead of the release (a development
    % copy) is not an update.
    needupdate = pkg.i_isnewerversion(v2web, v1local);
catch ME
    if needed(3) && isempty(v2web)
        % The released version was never retrieved, so "up to date" cannot
        % be asserted. Let the caller report the failure.
        rethrow(ME);
    end
    disp(ME.message);
end
    % if nargout > 1 && needed(2), v1local = v1local; end

if nargout > 3 && needed(4)
        im = [];
    end
end
