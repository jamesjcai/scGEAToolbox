function callback_CheckUpdates(~, ~)

% Col = webread('https://api.github.com/repos/jamesjcai/scGEAToolbox')
[needupdate, v1, v2] = pkg.i_vercheck;
% Check for updates.
if needupdate
    answer = questdlg(sprintf('There is a new version of scGEAToolbox (%s vs. %s). Learn how to install?', v2, v1));
    if strcmp(answer, 'Yes')
        web('https://scgeatoolbox.readthedocs.io/en/latest/quick_installation.html');
    end
else
    helpdlg(sprintf('scGEAToolbox (%s) is up to date.', v1), '');
end
end