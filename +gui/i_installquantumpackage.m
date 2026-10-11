function i_installquantumpackage(src, ~)
% I_INSTALLQUANTUMPACKAGE  Install or upgrade the MATLAB Quantum Computing Support Package.
%   Reports the installed version and offers to open the package's page in
%   Add-On Explorer, which shows Install or Update as appropriate. There is
%   no public feed of the latest version: mathworks/Quantum-Computing-MATLAB
%   publishes no GitHub releases, and the File Exchange version history
%   lags the actual release, so Add-On Explorer is the only reliable judge.

[parentfig, ~] = gui.gui_getfigsce(src);

ADDON_NAME = 'MATLAB Support Package for Quantum Computing';
% Support-package base code, from the add-on's .mlpkginstall signpost.
% Add-On Explorer opens a support package's page by this code, not by
% its File Exchange uuid.
ADDON_ID = 'QUANTUM';
FEX_URL = 'https://www.mathworks.com/matlabcentral/fileexchange/125425-matlab-support-package-for-quantum-computing';

installedVersion = i_getInstalledVersion(ADDON_ID);
if isempty(installedVersion)
    msg = sprintf(['%s is not installed.\n' ...
        'Open it in Add-On Explorer to install?'], ADDON_NAME);
    dlgTitle = 'Install Package';
else
    msg = sprintf(['%s v%s is installed.\n' ...
        'Open it in Add-On Explorer to check for an update?'], ...
        ADDON_NAME, installedVersion);
    dlgTitle = 'Package Update';
end
answer = gui.myQuestdlg(parentfig, msg, dlgTitle, ...
    {'Add-On Explorer', 'Cancel'}, 'Add-On Explorer');
if strcmp(answer, 'Add-On Explorer')
    try
        %#exclude matlab.internal.addons.launchers.showExplorer
        matlab.internal.addons.launchers.showExplorer('scgeatool', ...
            'identifier', ADDON_ID);
    catch
        % Internal launcher unavailable (e.g. deployed); use the web page.
        web(FEX_URL, '-browser');
    end
end
end


% -----------------------------------------------------------------------------
function version = i_getInstalledVersion(identifier)
% Return the installed version, or '' when the add-on is not installed.
version = '';
if ismcc || isdeployed
    return;
end
try
    addons = matlab.addons.installedAddons;
    idx = find(addons.Identifier == identifier, 1);
    if ~isempty(idx)
        version = char(addons.Version(idx));
    end
catch
    % Add-on query unavailable; treat as not installed.
end
end
