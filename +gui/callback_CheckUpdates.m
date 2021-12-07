function callback_CheckUpdates(~,~)
[needupdate]=pkg.i_vercheck;
% Check for updates.
if needupdate
    answer = questdlg('There is a new version of scGEAToolbox. Learn how to install?');
    if strcmp(answer,'Yes')
        web('https://scgeatoolbox.readthedocs.io/en/latest/quick_installation.html');
    end
else
    helpdlg('scGEAToolbox is up to date.','')
end
end