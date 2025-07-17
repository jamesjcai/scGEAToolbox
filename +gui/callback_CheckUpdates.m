function callback_CheckUpdates(src, ~)    

    toolboxPath = fileparts(fileparts(mfilename('fullpath')));
    %versionfile = fullfile(toolboxPath, 'VERSION.mat');
    %if exist(versionfile,'file')
    %    delete(versionfile);
    %end
    [FigureHandle] = gui.gui_getfigsce(src);
    % Col = webread('https://api.github.com/repos/jamesjcai/scGEAToolbox')
    try
        [majneedupdate, v_old, v_new] = pkg.i_majvercheck;
    catch ME
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier)
        % 'Could not access server.'
        return;
    end

    if majneedupdate
        if ~pkg.e_runningasaddons || ismcc || isdeployed || isa(src, 'matlab.apps.AppBase')
            answer = gui.myQuestdlg(FigureHandle, ...
                sprintf(['There is a new version of scGEAToolbox ' ...
                '(%s vs. %s). Learn how to upgrade?'], ...
                v_new, v_old));
            if strcmp(answer, 'Yes')
                gui.gui_showrefinfo('Quick Installation', FigureHandle);
            end
        else
            answer = gui.myQuestdlg(FigureHandle, ...
                sprintf('There is a new version of scGEAToolbox (%s vs. %s). Upgrade?', ...
                v_new, v_old));
            if strcmp(answer, 'Yes')
                try
                    fw = gui.myWaitbar(FigureHandle);
                    gui.myWaitbar(FigureHandle, fw, false, '', 'Downloading...', 1/5);
                    instURL = 'https://api.github.com/repos/jamesjcai/scGEAToolbox/releases/latest';
                    %[~, instName] = fileparts(fileparts(fileparts(instURL)));
                    instRes = webread(instURL);
                    % fprintf('Downloading %s %s\n', instName, instRes.name);
                    % websave(instRes.assets.name, instRes.assets.browser_download_url);
    
                    toolboxURL = instRes.assets.browser_download_url;
                    tempZip = fullfile(tempdir, instRes.assets.name);
    
                    %toolboxURL = sprintf('https://github.com/jamesjcai/scGEAToolbox/releases/download/v%s/scGEAToolbox.mltbx', v2);
                    %tempZip = fullfile(tempdir, "ToolboxUpdate.mltbx");
                    websave(tempZip, toolboxURL);
                    gui.myWaitbar(FigureHandle, fw, false, '', 'Installing...', 2/5);
                    warning off
                    matlab.addons.install(tempZip);
                    gui.myWaitbar(FigureHandle, fw, false, '', 'Post-install processing...', 3/5);
                    pause(2)
    
                    versionfile = fullfile(toolboxPath, 'VERSION.mat');
                    if exist(versionfile,'file')
                        delete(versionfile);
                    end
                    
                catch ME
                    gui.myErrordlg(FigureHandle, ME.message,'');
                    return;
                end
                warning on
                gui.myWaitbar(FigureHandle, fw, false, '', 'Update complete!', 5/5);
                pause(3)
                gui.myWaitbar(FigureHandle, fw);
                if strcmp('Yes', gui.myQuestdlg(FigureHandle, 'Restart SCGEATOOL?'))
                    close(FigureHandle);
                    pause(1);
                    run("scgeatool.m");
                end
            end
        end
    else
        gui.myHelpdlg(FigureHandle, sprintf('scGEAToolbox (%s) is up to date.', v_old));
        %answer=gui.myQuestdlg(FigureHandle, 'Check for minor updates?','');
        %if strcmp(answer,'Yes')
        %    pkg.i_minvercheck(FigureHandle);
        %end
    end

end