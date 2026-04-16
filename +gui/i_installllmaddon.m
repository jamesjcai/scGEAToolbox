function [done] = i_installllmaddon(src, ~)
% I_INSTALLLLMADDON  Install or upgrade the LLMs-with-MATLAB File Exchange addon.
%   Checks the installed version against the latest GitHub release. If an
%   upgrade is available (or the addon is absent), offers to auto-install the
%   .mltbx from the GitHub release, or falls back to opening the File Exchange
%   page.

[parentfig, ~] = gui.gui_getfigsce(src);
done = false;

ADDON_NAME = 'Large Language Models (LLMs) with MATLAB';
GITHUB_REPO = 'matlab-deep-learning/llms-with-matlab';
FEX_URL = ['https://www.mathworks.com/matlabcentral/fileexchange/' ...
           '163796-large-language-models-llms-with-matlab'];

% --- Installed state ---------------------------------------------------------
[isInstalled, installedVersion] = i_getInstalledVersion(ADDON_NAME);

% --- Latest version from GitHub ----------------------------------------------
[latestVersion, mltbxUrl] = i_fetchLatestRelease(GITHUB_REPO);

% --- Decide what to do -------------------------------------------------------
if isInstalled
    if ~isempty(latestVersion) && ~strcmp(installedVersion, latestVersion)
        answer = gui.myQuestdlg(parentfig, ...
            sprintf(['%s is installed (v%s).\n' ...
                     'A newer version (v%s) is available. Upgrade now?'], ...
                     ADDON_NAME, installedVersion, latestVersion), ...
            'Addon Update', {'Upgrade', 'Open File Exchange', 'Cancel'}, 'Upgrade');
    else
        if isempty(latestVersion)
            msg = sprintf(['%s (v%s) is installed.\n' ...
                '(Could not check for updates — network unavailable.)'], ...
                ADDON_NAME, installedVersion);
        else
            msg = sprintf('%s (v%s) is up to date.', ADDON_NAME, installedVersion);
        end
        gui.myHelpdlg(parentfig, msg);
        done = true;
        return;
    end
else
    if ~isempty(latestVersion)
        installMsg = sprintf('%s is not installed. Install v%s now?', ...
            ADDON_NAME, latestVersion);
    else
        installMsg = sprintf('%s is not installed. Install now?', ADDON_NAME);
    end
    answer = gui.myQuestdlg(parentfig, installMsg, ...
        'Install Addon', {'Install', 'Open File Exchange', 'Cancel'}, 'Install');
end

switch answer
    case {'Install', 'Upgrade'}
        done = i_doInstall(parentfig, mltbxUrl, FEX_URL);
    case 'Open File Exchange'
        web(FEX_URL, '-browser');
end
end


% -----------------------------------------------------------------------------
function [isInstalled, version] = i_getInstalledVersion(addonName)
isInstalled = false;
version = '';
try
    pkgs = matlab.addons.installedAddons;
    if isempty(pkgs), return; end
    idx = strcmp(pkgs.Name, addonName);
    if any(idx)
        isInstalled = true;
        v = pkgs.Version(idx);
        if iscell(v)
            version = char(v{1});
        else
            version = char(v);
        end
    end
catch
    % Addon query failed — treat as not installed
end
end


% -----------------------------------------------------------------------------
function [version, mltbxUrl] = i_fetchLatestRelease(githubRepo)
% Query the GitHub releases API for the latest tag and any .mltbx asset.
version = '';
mltbxUrl = '';
try
    apiUrl = sprintf('https://api.github.com/repos/%s/releases/latest', githubRepo);
    opts = weboptions('Timeout', 10, 'ContentType', 'json');
    data = webread(apiUrl, opts);

    tag = data.tag_name;
    if startsWith(tag, 'v')
        tag = tag(2:end);
    end
    version = tag;

    % Look for a .mltbx asset in the release
    assets = data.assets;
    if isstruct(assets) && ~isempty(assets)
        for k = 1:numel(assets)
            name = string(assets(k).name);
            if endsWith(name, '.mltbx')
                mltbxUrl = char(assets(k).browser_download_url);
                break;
            end
        end
    end
catch
    % Network unavailable, no releases, or unexpected JSON shape
end
end


% -----------------------------------------------------------------------------
function done = i_doInstall(parentfig, mltbxUrl, fexUrl)
done = false;
if ~isempty(mltbxUrl)
    fw = gui.myWaitbar(parentfig);
    try
        tempFile = fullfile(tempdir, 'llms_with_matlab.mltbx');
        websave(tempFile, mltbxUrl);
        gui.myWaitbar(parentfig, fw);
        matlab.addons.install(tempFile, true, 'overwrite');
        done = true;
        gui.myHelpdlg(parentfig, ...
            ['Installation complete. ' ...
             'You may need to restart MATLAB for changes to take effect.']);
    catch ME
        gui.myWaitbar(parentfig, fw);
        gui.myErrordlg(parentfig, ['Installation failed: ' ME.message]);
        web(fexUrl, '-browser');
    end
else
    % No .mltbx asset in the release — open the File Exchange page
    gui.myHelpdlg(parentfig, ...
        'No direct download found. Opening the File Exchange page to install manually.');
    web(fexUrl, '-browser');
end
end
