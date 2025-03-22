function callback_InstallApp(src, ~)

[FigureHandle] = gui.gui_getfigsce(src);
toolboxPath = fileparts(fileparts(mfilename('fullpath')));
filein = fullfile(toolboxPath, 'scgeatoolApp.mlappinstall');

if exist(filein, "file")
    try
        appinfo = matlab.apputil.install(filein);
        % assignin("base", "a", appinfo);
        % Define the app name and usage instructions
        appName = appinfo.name;
        appUsage = sprintf(['To open the app, follow these steps:\n\n',...
            '1. Go to the "Apps" tab in the MATLAB toolstrip.\n',...
            '2. Look for "%s" in the list of installed apps.\n',...
            '3. Click on it to launch.'], appName);
        gui.myHelpdlg(FigureHandle, appUsage, 'Installation Successful');
    catch ME
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    end
else
    gui.myHelpdlg(FigureHandle, 'Installation file is missing.', '');
end