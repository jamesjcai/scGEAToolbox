function [done] = i_setwrkdir(preftagname)
%I_SETWRKDIR - set workding directory
%see also: I_SETPYENV, I_SETRENV 
[done] = false;

% preftagname = 'externalwrkpath';    i_setextwd
% preftagname = 'netanalywrkpath';    i_setnetwd

if ~ispref('scgeatoolbox', preftagname)
    answer = questdlg('Working directory has not been set up. Locate a folder?');
    if ~strcmp(answer, 'Yes'), return; end
    if ispc
        [~,b]=system("echo %username%");
        pathdefult = sprintf('C:\\Users\\%s\\Documents\\', ...
            string(deblank(b)));
    else
        pathdefult = '';
    end
    [done] = ix_setwdpath(pathdefult);
else
    s = getpref('scgeatoolbox', preftagname);
    answer1 = questdlg(sprintf('%s', s), ...
        'Working Root', ...
        'Use this', 'Use another', 'Cancel', 'Use this');
    switch answer1
        case 'Use this'
            done = true;
        case 'Use another'
            [done] = ix_setwdpath(s);
    end
end



    function [done] = ix_setwdpath(deflt)
        done = false;
        answer=questdlg('Where to save working files?','','Use Temporary Folder', ...
            'Select a Folder','Cancel','Use Temporary Folder');
        switch answer
            case 'Select a Folder'
                [seltpath] = uigetdir(deflt);
                if seltpath==0, return; end
                if ~isfolder(seltpath), return; end
            case 'Use Temporary Folder'
                seltpath = tempdir;
            case 'Cancel'
                return;
        end
        disp(['User selected: ', seltpath]);
        try
            setpref('scgeatoolbox', preftagname, seltpath);
        catch ME
            errordlg(ME.message);
            return;
        end
        done = true;
    end

end
