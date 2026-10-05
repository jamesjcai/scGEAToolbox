function [done] = i_installllmaddon(src, ~)
% I_INSTALLLLMADDON  Install or upgrade the LLMs-with-MATLAB add-on.
%   Installs the File Exchange package llms_with_matlab with the MATLAB
%   package manager (mpminstall/mpmupdate), which resolves the latest
%   version itself. That replaced downloading a .mltbx from the GitHub
%   release, which stopped working when the releases stopped carrying one
%   (v4.9.0 has none). If the add-on was installed the old way, as a
%   .mltbx through the Add-On Explorer, the upgrade uninstalls that copy
%   first so the path does not hold two. Falls back to the File Exchange
%   page when the package manager cannot do it.

[parentfig, ~] = gui.gui_getfigsce(src);
done = false;

PKG_NAME = "llms_with_matlab";
ADDON_NAME = "Large Language Models (LLMs) with MATLAB";
FEX_URL = "https://www.mathworks.com/matlabcentral/fileexchange/" + ...
    "163796-large-language-models-llms-with-matlab";
LICENSE_NOTE = "The add-on comes under its own license terms, " + ...
    "which are placed in its installation folder.";

% --- Installed state and latest version -------------------------------------
[installedVersion, isPackage, legacyId] = i_getInstalled(PKG_NAME, ADDON_NAME);
latestVersion = i_getLatestVersion(PKG_NAME, isPackage);

% --- Decide what to do -------------------------------------------------------
if installedVersion ~= ""
    if latestVersion ~= "" && latestVersion ~= installedVersion
        answer = gui.myQuestdlg(parentfig, ...
            sprintf(['%s is installed (v%s).\n' ...
                     'A newer version (v%s) is available. Upgrade now?\n\n%s'], ...
                     ADDON_NAME, installedVersion, latestVersion, LICENSE_NOTE), ...
            'Addon Update', {'Upgrade', 'Open File Exchange', 'Cancel'}, 'Upgrade');
    else
        if latestVersion == ""
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
    if latestVersion ~= ""
        installMsg = sprintf('%s is not installed. Install v%s now?\n\n%s', ...
            ADDON_NAME, latestVersion, LICENSE_NOTE);
    else
        installMsg = sprintf('%s is not installed. Install now?\n\n%s', ...
            ADDON_NAME, LICENSE_NOTE);
    end
    answer = gui.myQuestdlg(parentfig, installMsg, ...
        'Install Addon', {'Install', 'Open File Exchange', 'Cancel'}, 'Install');
end

switch answer
    case {'Install', 'Upgrade'}
        done = i_doInstall(parentfig, PKG_NAME, isPackage, legacyId, FEX_URL);
    case 'Open File Exchange'
        web(FEX_URL, '-browser');
    otherwise
        % Cancel, or the dialog was closed: nothing to do.
end
end


% -----------------------------------------------------------------------------
function [version, isPackage, legacyId] = i_getInstalled(pkgName, addonName)
% version is "" when the add-on is absent. isPackage is true for a
% package-manager install; legacyId is the Add-On identifier of a .mltbx
% install, which the package manager does not list, and "" otherwise.
version = "";
isPackage = false;
legacyId = "";
try
    pkg = mpmlist(Name=pkgName);
    if ~isempty(pkg)
        version = string(pkg(1).Version);
        isPackage = true;
        return;
    end
    addons = matlab.addons.installedAddons;
    idx = find(addons.Name == addonName, 1);
    if ~isempty(idx)
        version = string(addons.Version(idx));
        legacyId = string(addons.Identifier(idx));
    end
catch
    % Query failed: treat as not installed
end
end


% -----------------------------------------------------------------------------
function version = i_getLatestVersion(pkgName, isPackage)
% Ask the package manager what it would install, without installing it.
% "" when the repository cannot be reached.
version = "";
try
    if isPackage
        pkg = mpmupdate(pkgName, DryRun=true, Prompt=false, Verbosity="quiet");
    else
        pkg = mpminstall(pkgName, DryRun=true, Prompt=false, Verbosity="quiet");
    end
    version = string(pkg(1).Version);
catch
    % Network unavailable or the package is not in the repository
end
end


% -----------------------------------------------------------------------------
function done = i_doInstall(parentfig, pkgName, isPackage, legacyId, fexUrl)
done = false;
removedLegacy = false;
fw = gui.myWaitbar(parentfig);
try
    if isPackage
        pkg = mpmupdate(pkgName, Prompt=false, Verbosity="quiet");
    else
        % The package is installed first and the old .mltbx copy removed
        % only if it gets in the way. It used to be removed first, so a
        % failed install left neither.
        try
            pkg = mpminstall(pkgName, Prompt=false, Verbosity="quiet");
        catch firstErr
            if legacyId == "", rethrow(firstErr); end
            matlab.addons.uninstall(legacyId);
            removedLegacy = true;
            pkg = mpminstall(pkgName, Prompt=false, Verbosity="quiet");
        end
    end
    gui.myWaitbar(parentfig, fw);
    done = true;
    gui.myHelpdlg(parentfig, sprintf('Installed %s v%s.', ...
        pkgName, string(pkg(1).Version)));
catch ME
    gui.myWaitbar(parentfig, fw);
    msg = ['Installation failed: ' ME.message];
    if removedLegacy
        msg = [msg newline 'The previously installed copy was removed ' ...
            'to make way for it; reinstall it from the page that opens.'];
    end
    gui.myErrordlg(parentfig, msg);
    web(fexUrl, '-browser');
end
end
