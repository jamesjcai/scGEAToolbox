function callback_CheckUpdates(src, ~)
% Check for updates.
FigureHandle = src.Parent.Parent;
% Col = webread('https://api.github.com/repos/jamesjcai/scGEAToolbox')
[majneedupdate, v1, v2] = pkg.i_majvercheck;

if majneedupdate
    answer = questdlg(sprintf('There is a new version of scGEAToolbox (%s vs. %s). Learn how to install?', v2, v1));
    if strcmp(answer, 'Yes')
        gui.gui_uishowrefinfo('Quick Installation', FigureHandle);

        %web('https://scgeatoolbox.readthedocs.io/en/latest/quick_installation.html');        
        % prompt = {'Copy the following code and run it in MATLAB:'};
        % dlgtitle = 'Quick Installation';
        % fieldsize = [18 75];
        % definput = {sprintf('tic;\ndisp(''Installing scGEAToolbox...'')\nunzip(''https://github.com/jamesjcai/scGEAToolbox/archive/main.zip'');\naddpath(''./scGEAToolbox-main'');\ntoc;\nif exist(''scgeatool.m'',''file'')\n    disp(''scGEAToolbox installed!'')\nend\nsavepath(fullfile(userpath,''pathdef.m''));')};
        % inputdlg(prompt,dlgtitle,fieldsize,definput);
    end
else
    waitfor(helpdlg(sprintf('scGEAToolbox (%s) is up to date.', v1), ''));
    answer=questdlg('Check for minor updates?','');
    if strcmp(answer,'Yes')
        pkg.i_minvercheck(FigureHandle);
    end
end
