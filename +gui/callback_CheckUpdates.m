function callback_CheckUpdates(src, ~)    

    toolboxPath = fileparts(fileparts(mfilename('fullpath')));
    %versionfile = fullfile(toolboxPath, 'VERSION.mat');
    %if exist(versionfile,'file')
    %    delete(versionfile);
    %end
    [FigureHandle] = gui.gui_getfigsce(src);
    % Col = webread('https://api.github.com/repos/jamesjcai/scGEAToolbox')
    try
        [majneedupdate, v_old, v_new] = pkg.i_majvercheck(true);
    catch ME
        gui.myErrordlg(FigureHandle, 'Could not access server.');
        return;
    end
    if majneedupdate
        if ~pkg.e_runningasaddons || ismcc || isdeployed
            answer = gui.myQuestdlg(FigureHandle, sprintf('There is a new version of scGEAToolbox (%s vs. %s). Learn how to upgrade?', ...
                v_new, v_old));
            if strcmp(answer, 'Yes')               
                gui.gui_showrefinfo('Quick Installation', FigureHandle);
            end
        else
            answer = gui.myQuestdlg(FigureHandle, sprintf('There is a new version of scGEAToolbox (%s vs. %s). Upgrade?', ...
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
            %{
            toolboxPath = fileparts(fileparts(mfilename('fullpath')));
            if isempty(toolboxPath)
                gui.myErrordlg(FigureHandle, "Toolbox not found on MATLAB path!", "Update Error");
                return;
            end

            if isfile(fullfile(toolboxPath, 'DEVMODE.txt'))
                gui.myErrordlg(FigureHandle, "Upgrade is disabled in development mode.", "Update Blocked");
                return;
            end

            if upgradeToolbox
                if strcmp('Yes', gui.myQuestdlg(FigureHandle, 'Restart SCGEATOOL?'))
                    close(FigureHandle);
                    pause(1);
                    run("scgeatool.m");
                end
            end
            %}
        end
    else
        gui.myHelpdlg(FigureHandle, sprintf('scGEAToolbox (%s) is up to date.', v_old));
        %answer=gui.myQuestdlg(FigureHandle, 'Check for minor updates?','');
        %if strcmp(answer,'Yes')
        %    pkg.i_minvercheck(FigureHandle);
        %end
    end

    %{
    function [done] = upgradeToolbox()
        done = false;
        backupDir = fullfile(tempdir, "ToolboxBackup"); % Backup location
        % zipURL = sprintf('https://github.com/jamesjcai/scGEAToolbox/releases/download/v%s/scGEAToolbox.mltbx', v2); % Update with actual URL
        zipURL = sprintf('https://github.com/jamesjcai/scGEAToolbox/archive/refs/tags/v%s.zip',v2);
        tempZip = fullfile(tempdir, "ToolboxUpdate.zip");
    
%        zipURL = sprintf('https://github.com/jamesjcai/scGEAToolbox/releases/download/v%s/scGEAToolbox.mltbx', v2); % Update with actual URL
%        tempZip = fullfile(tempdir, "ToolboxUpdate.mltbx");
%        websave(tempZip, zipURL);    
%        matlab.addons.install(tempZip);


        try
            % Step 1: Backup the existing toolbox
            disp("Backing up current version...");
            if exist(backupDir, 'dir')
                rmdir(backupDir, 's'); % Remove old backup
            end
            copyfile(toolboxPath, backupDir); % Copy toolbox to backup location
            
            % Step 2: Download latest toolbox ZIP
            disp("Downloading update...");
            websave(tempZip, zipURL);
    
            % Step 3: Extract ZIP to toolbox location (overwrite existing files)
            disp("Extracting files...");
            unzip(tempZip, toolboxPath);
            
            % Step 4: Clean up temporary files
            delete(tempZip);
    
            % Step 5: Refresh MATLAB path
            rehash toolboxcache;
            savepath;

            versionfile = fullfile(toolboxPath, 'VERSION.mat');
            if exist(versionfile,'file')
                delete(versionfile);
            end
            % Step 6: Notify user
            disp("Toolbox updated successfully!");
            msgbox("Update complete! Restart MATLAB to apply changes.", "Update Successful");
            done = true;
        catch ME
            % Step 7: Rollback on failure
            warning(ME.identifier, 'Update failed: %s', ME.message);
            msgbox("Update failed. Restoring previous version...", "Update Error", "error");
    
            % Restore backup
            if exist(backupDir, 'dir')
                disp("Restoring previous version...");
                rmdir(toolboxPath, 's'); % Delete broken update
                copyfile(backupDir, toolboxPath); % Restore backup
                msgbox("Rollback complete! Restart MATLAB to apply changes.", "Rollback Successful");
            else
                msgbox("Rollback failed: No backup found.", "Rollback Error", "error");
            end
        end
    end
    %}
    
end