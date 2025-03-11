function [done] = i_setextwd(src, ~)
%I_SETEXTWD - set externel workding directory

%see also: I_SETPYENV, I_SETRENV 
[FigureHandle, ~] = gui.gui_getfigsce(src);
[done] = gui.i_setwrkdir('externalwrkpath');
if done
     gui.myHelpdlg(FigureHandle, "External program working directory is set successfully.", '');
end


%{
[done] = false;

preftagname = 'externalwrkpath';

if ~ispref('scgeatoolbox', preftagname)
    answer = gui.myQuestdlg(FigureHandle, 'External working directory has not been set up. Locate a folder?');
    if ~strcmp(answer, 'Yes'), return; end
    if ispc
        [~,b]=system("echo %username%");
        pathdefult = sprintf('C:\Users\%s\Documents',string(deblank(b)));
    else
        pathdefult = '';
    end
    [done] = ix_setrenv(pathdefult);
    % setpref('scgeatoolbox','rexecutablepath',s);
else
    s = getpref('scgeatoolbox', preftagname);
    answer = gui.myQuestdlg(FigureHandle, sprintf('%s', s), ...
        'External Working Directory', ...
        'Use this', 'Use another', 'Cancel', 'Use this');
    switch answer
        case 'Use this'
            done = true;
        case 'Use another'
            [done] = ix_setrenv(s);
    end
end

if done
     gui.myHelpdlg(FigureHandle, "External working directory is set successfully.", '');
end


    function [done] = ix_setrenv(deflt)
        % selpath = uigetdir;
        done = false;
        [seltpath] = uigetdir(deflt);
        if isequal(seltpath, 0)
            %disp('User selected Cancel');
            return;
        else
            disp(['User selected: ', seltpath]);
            % fullfile(path)
            try
                setpref('scgeatoolbox', preftagname, seltpath);
            catch ME
                errordlg(ME.message);
                return;
            end
            done = true;
        end
    end

end
%}