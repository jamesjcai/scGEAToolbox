function i_installjava(src, ~)
% I_INSTALLJAVA  Install and configure an OpenJDK runtime for MATLAB.
%   From R2026b MATLAB no longer ships a JRE, so features that need Java
%   fail until one is configured with jenv. This checks the current Java
%   environment and, if none is set up, either opens the "MATLAB Support
%   for OpenJDK" add-on in Add-On Explorer or downloads an Eclipse Adoptium
%   JRE into the toolbox add-on folder and points jenv at it. The new
%   runtime is picked up after MATLAB restarts.
%
%   See https://www.mathworks.com/help/matlab/matlab_external/configure-your-system-to-use-java.html

[parentfig, ~] = gui.gui_getfigsce(src);

JAVA_VERSION = 21;
ADOPTIUM_URL = 'https://adoptium.net/';
% Support-package base code, from the add-on's .mlpkginstall signpost.
% Add-On Explorer opens a support package's page by this code, not by
% its File Exchange uuid.
ADDON_ID = 'MLJRE';
ADDON_URL = 'https://www.mathworks.com/matlabcentral/fileexchange/180823-matlab-support-for-openjdk';

% jenv is a built-in (R2021b+), so exist(..., 'file') cannot detect it.
if isMATLABReleaseOlderThan("R2021b")
    gui.myHelpdlg(parentfig, ['This MATLAB release has no jenv function. ' ...
        'It uses the Java runtime shipped with MATLAB.']);
    return;
end

% --- Current state -----------------------------------------------------------
je = jenv;
if je.Status == "loaded"
    gui.myHelpdlg(parentfig, sprintf('Java is loaded:\n%s\n\nHome: %s', ...
        je.Version, je.Home));
    return;
end
if strlength(je.Configuration) > 0 && je.Configuration ~= "system" ...
        && isfolder(je.Configuration)
    gui.myHelpdlg(parentfig, sprintf(['MATLAB is configured to use the Java at\n%s\n' ...
        'but it is not loaded yet. Restart MATLAB to load it.'], je.Configuration));
    return;
end
if i_isaddoninstalled(ADDON_ID)
    gui.myHelpdlg(parentfig, ['MATLAB Support for OpenJDK is installed but ' ...
        'Java is not loaded. Run jenv -clear, then restart MATLAB.']);
    return;
end

% --- Offer to install ----------------------------------------------------------
downloadLabel = sprintf('Download OpenJDK %d', JAVA_VERSION);
answer = gui.myQuestdlg(parentfig, ...
    sprintf(['Java is not configured for MATLAB. Features that call Java ' ...
    'will not work until it is.\n\n' ...
    'Add-On Explorer: install the "MATLAB Support for OpenJDK" add-on ' ...
    '(MathWorks-supported; installs OpenJDK 8 on Windows/Linux, 11 on Mac).\n\n' ...
    '%s: fetch the Eclipse Adoptium JRE directly and configure jenv ' ...
    '(no MathWorks sign-in needed).\n\nEither way, restart MATLAB afterwards.'], ...
    downloadLabel), ...
    'Install Java', {'Add-On Explorer', downloadLabel, 'Cancel'}, 'Add-On Explorer');

switch answer
    case 'Add-On Explorer'
        try
            %#exclude matlab.internal.addons.launchers.showExplorer
            matlab.internal.addons.launchers.showExplorer('scgeatool', ...
                'identifier', ADDON_ID);
        catch
            % Internal launcher unavailable (e.g. deployed); use the web page.
            web(ADDON_URL, '-browser');
        end
        gui.myHelpdlg(parentfig, ['Click Install in the Add-On Explorer, ' ...
            'then run jenv -clear and restart MATLAB.']);
    case downloadLabel
        fw = gui.myWaitbar(parentfig);
        try
            javaHome = i_downloadjre(JAVA_VERSION);
            jenv(javaHome);
            gui.myWaitbar(parentfig, fw);
            gui.myHelpdlg(parentfig, sprintf(['OpenJDK %d installed in\n%s\n\n' ...
                'Restart MATLAB to load it.'], JAVA_VERSION, javaHome));
        catch ME
            gui.myWaitbar(parentfig, fw);
            gui.myErrordlg(parentfig, sprintf(['Java installation failed: %s\n\n' ...
                'Install OpenJDK from %s, then run jenv("<JRE folder>").'], ...
                ME.message, ADOPTIUM_URL));
        end
    otherwise
        % Cancelled or dialog closed: nothing to do.
end
end


% -----------------------------------------------------------------------------
function javaHome = i_downloadjre(javaVersion)
% Download and extract an Adoptium JRE; return the folder that holds bin/java.
switch computer('arch')
    case 'win64'
        os = 'windows';
        arch = 'x64';
        ext = '.zip';
    case 'maci64'
        os = 'mac';
        arch = 'x64';
        ext = '.tar.gz';
    case 'maca64'
        os = 'mac';
        arch = 'aarch64';
        ext = '.tar.gz';
    case 'glnxa64'
        os = 'linux';
        arch = 'x64';
        ext = '.tar.gz';
    otherwise
        error('No OpenJDK download is defined for platform %s.', computer('arch'));
end

if strcmp(os, 'linux') && javaVersion == 21
    % MATLAB on Linux supports OpenJDK 21 only up to 21.0.2+13.
    url = ['https://github.com/adoptium/temurin21-binaries/releases/download/' ...
        'jdk-21.0.2%2B13/OpenJDK21U-jre_x64_linux_hotspot_21.0.2_13.tar.gz'];
else
    url = sprintf(['https://api.adoptium.net/v3/binary/latest/%d/ga/' ...
        '%s/%s/jre/hotspot/normal/eclipse'], javaVersion, os, arch);
end

installBase = fullfile(prefdir, 'scgeatoolbox_addons', sprintf('openjdk%d', javaVersion));
if isfolder(installBase)
    rmdir(installBase, 's');
end
mkdir(installBase);

archiveFile = fullfile(tempdir, ['openjdk_jre' ext]);
cleanupArchive = onCleanup(@() i_deletefile(archiveFile));
websave(archiveFile, url, weboptions('Timeout', 120));
if strcmp(ext, '.zip')
    unzip(archiveFile, installBase);
else
    untar(archiveFile, installBase);
end

javaHome = i_findjavahome(installBase);
if isempty(javaHome)
    error('The downloaded archive does not contain a Java runtime.');
end
end


% -----------------------------------------------------------------------------
function javaHome = i_findjavahome(folder)
% Return the folder whose bin/ holds the java executable (Contents/Home on macOS).
javaHome = '';
if ispc
    exeName = 'java.exe';
else
    exeName = 'java';
end
hits = dir(fullfile(folder, '**', 'bin', exeName));
if ~isempty(hits)
    javaHome = fileparts(hits(1).folder);
end
end


% -----------------------------------------------------------------------------
function tf = i_isaddoninstalled(identifier)
tf = false;
try
    addons = matlab.addons.installedAddons;
    tf = any(addons.Identifier == identifier);
catch
    % Add-on query unavailable (e.g. deployed); treat as not installed.
end
end


% -----------------------------------------------------------------------------
function i_deletefile(fileName)
if isfile(fileName)
    delete(fileName);
end
end
