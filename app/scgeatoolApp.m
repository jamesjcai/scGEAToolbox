function scgeatoolApp

disp('scgeatoolApp is running.')
a = which('scgeatool');
if ~isempty(a)
    y0 = true;
else
    y0 = false;
end
y1 = contains(a, 'MATLAB Add-Ons');

if y0 && ~y1
    helpdlg(['scGEATool has been installed as development package ' ...
        '(instead of an Add-On). Run scgeatool manually.'],'');
    return;
end
t = matlab.addons.installedAddons;
[y2,b] = ismember('scGEAToolbox', t.Name);
if y2
    y3 = t.Enabled(b);
else
    y3 = false;
end
if y0 && y1 && y2 && y3
    scgeatool
    disp('scgeatool starts. Enjoy exploring!')
end


if ~y0 && ~y1 && ~y2
    if strcmp('Yes', questdlg('Install scGEAToolbox Add-on?',''))
        try
            instURL = 'https://api.github.com/repos/jamesjcai/scGEAToolbox/releases/latest';
            [~, instName] = fileparts(fileparts(fileparts(instURL)));
            instRes = webread(instURL);
            fprintf('Downloading %s %s ...... ', instName, instRes.tag_name);
           
            toolboxURL = instRes.assets.browser_download_url;
            tempZip = fullfile(tempdir, instRes.assets.name);
            websave(tempZip, toolboxURL);
            fprintf('Done.\n');
            
            fprintf('Installing ......');
            warning off
            matlab.addons.install(tempZip);
            fprintf('Done.\n');
        catch ME
            errordlg(ME.message, ME.identifier);
            return;
        end

        if strcmp('Yes', questdlg('Start scgeatool?',''))
            scgeatool;
        end
    end
end
