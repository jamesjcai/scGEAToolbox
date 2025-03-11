function callback_CheckUpdates(src, ~)    

    toolboxPath = fileparts(fileparts(mfilename('fullpath')));
    [FigureHandle] = gui.gui_getfigsce(src);
    % Col = webread('https://api.github.com/repos/jamesjcai/scGEAToolbox')
    try
        [majneedupdate, v_old, v_new] = pkg.i_majvercheck;
    catch ME
        errordlg('Could not access server.');
        return;
    end
    if majneedupdate
        if ~pkg.e_runningasaddons
            answer = gui.myQuestdlg(FigureHandle, sprintf('There is a new version of scGEAToolbox (%s vs. %s). Learn how to upgrade?', ...
                v_new, v_old));
            if strcmp(answer, 'Yes')               
                gui.gui_uishowrefinfo('Quick Installation', FigureHandle);
            end
        else
            answer = gui.myQuestdlg(FigureHandle, sprintf('There is a new version of scGEAToolbox (%s vs. %s). Upgrade?', ...
                v_new, v_old));
            if strcmp(answer, 'Yes')
                try
                    gui.gui_waitbar_adv(fw,1/5,'Downloading...');
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
                    gui.gui_waitbar_adv(fw,2/5,'Installing...');
                    warning off
                    matlab.addons.install(tempZip);
                    gui.gui_waitbar_adv(fw,3/5,'Post-install processing...');
                    pause(2)
    
                    versionfile = fullfile(toolboxPath, 'VERSION.mat');
                    if exist(versionfile,'file')
                        delete(versionfile);
                    end
                    
                catch ME
                    errordlg(ME.message,'');
                    return;
                end
                warning on
                gui.gui_waitbar_adv(fw,5/5,'Update complete!');
                pause(3)
                gui.gui_waitbar_adv(fw);
                if strcmp('Yes', gui.myQuestdlg(FigureHandle, 'Restart SCGEATOOL?'))
                    close(FigureHandle);
                    pause(1);
                    run("scgeatool.m");
                end
            end
            %{
            toolboxPath = fileparts(fileparts(mfilename('fullpath')));
            if isempty(toolboxPath)
                errordlg("Toolbox not found on MATLAB path!", "Update Error");
                return;
            end

            if isfile(fullfile(toolboxPath, 'DEVMODE.txt'))
                errordlg("Upgrade is disabled in development mode.", "Update Blocked");
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