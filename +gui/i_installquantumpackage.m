function i_installquantumpackage(src, ~)
% I_INSTALLQUANTUMPACKAGE  Install or upgrade the MATLAB Quantum Computing Support Package.
%   Checks the installed version against the latest GitHub release and prompts
%   the user to upgrade if a newer version is available. Support packages
%   require the MathWorks installer, so the function opens the File Exchange
%   page rather than installing automatically.

[parentfig, ~] = gui.gui_getfigsce(src);

ADDON_NAME = 'MATLAB Support Package for Quantum Computing';
GITHUB_REPO = 'mathworks/Quantum-Computing-MATLAB';
FEX_URL = 'https://www.mathworks.com/matlabcentral/fileexchange/125425-matlab-support-package-for-quantum-computing';

% --- Installed state ---------------------------------------------------------
[isInstalled, installedVersion] = i_getInstalledVersion(ADDON_NAME);

% --- Latest version from GitHub ----------------------------------------------
latestVersion = i_fetchLatestRelease(GITHUB_REPO);

% --- Decide what to do -------------------------------------------------------
if isInstalled
    if ~isempty(latestVersion) && ~strcmp(installedVersion, latestVersion)
        answer = gui.myQuestdlg(parentfig, ...
            sprintf(['%s is installed (v%s).\n' ...
                     'A newer version (v%s) is available. Open the installer page?'], ...
                     ADDON_NAME, installedVersion, latestVersion), ...
            'Package Update', {'Open Installer Page', 'Cancel'}, 'Open Installer Page');
        if strcmp(answer, 'Open Installer Page')
            web(FEX_URL, '-browser');
        end
    else
        if isempty(latestVersion)
            msg = sprintf(['%s (v%s) is installed.\n' ...
                '(Could not check for updates — network unavailable.)'], ...
                ADDON_NAME, installedVersion);
        else
            msg = sprintf('%s (v%s) is up to date.', ADDON_NAME, installedVersion);
        end
        gui.myHelpdlg(parentfig, msg);
    end
else
    if ~isempty(latestVersion)
        installMsg = sprintf(['%s is not installed.\n' ...
            'Latest version: v%s. Open the installer page?'], ...
            ADDON_NAME, latestVersion);
    else
        installMsg = sprintf(['%s is not installed.\n' ...
            'Open the MATLAB Support Package page to download the installer?'], ...
            ADDON_NAME);
    end
    answer = gui.myQuestdlg(parentfig, installMsg, ...
        'Install Package', {'Open Installer Page', 'Cancel'}, 'Open Installer Page');
    if strcmp(answer, 'Open Installer Page')
        web(FEX_URL, '-browser');
    end
end
end


% -----------------------------------------------------------------------------
function [isInstalled, version] = i_getInstalledVersion(addonName)
isInstalled = false;
version = '';
if ismcc || isdeployed
    return;
end
try
    %#exclude matlabshared.supportpkg.getInstalled
    pkgs = matlabshared.supportpkg.getInstalled;
    if isempty(pkgs), return; end
    idx = strcmp({pkgs.Name}, addonName);
    if any(idx)
        isInstalled = true;
        v = pkgs(idx).Version;
        if iscell(v)
            version = char(v{1});
        else
            version = char(v);
        end
    end
catch
    % Support package query unavailable
end
end


% -----------------------------------------------------------------------------
function version = i_fetchLatestRelease(githubRepo)
% Query the GitHub releases API for the latest tag name.
version = '';
try
    apiUrl = sprintf('https://api.github.com/repos/%s/releases/latest', githubRepo);
    opts = weboptions('Timeout', 10, 'ContentType', 'json');
    data = webread(apiUrl, opts);
    tag = data.tag_name;
    if startsWith(tag, 'v')
        tag = tag(2:end);
    end
    version = tag;
catch
    % Network unavailable or no releases published
end
end
