function [done] = i_setwrkdir(preftagname, parentfig)
%I_SETWRKDIR - set workding directory
%see also: I_SETPYENV, I_SETRENV 

if nargin<2, parentfig = []; end
[done] = false;

% preftagname = 'externalwrkpath';    i_setextwd
% preftagname = 'netanalywrkpath';    i_setnetwd
if nargin < 1, preftagname = 'externalwrkpath'; end

if ~ispref('scgeatoolbox', preftagname)
    answer = gui.myQuestdlg(parentfig, ['Working directory has ' ...
        'not been set up. Locate a folder?']);
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
    answer1 = gui.myQuestdlg(parentfig, sprintf('%s', s), ...
        'Working Root', ...
        {'Use this', 'Use another', 'Cancel'}, 'Use this');
    switch answer1
        case 'Use this'
            done = true;
        case 'Use another'
            [done] = ix_setwdpath(s);
    end
end



    function [done] = ix_setwdpath(deflt)
        done = false;
        answer=gui.myQuestdlg(parentfig, 'Where to save working files?','',...
            {'Use Temporary Folder', ...
            'Select a Folder','Cancel'},'Use Temporary Folder');
        switch answer
            case 'Select a Folder'
                [seltpath] = uigetdir(deflt);
                if isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure')
                    figure(parentfig);
                end
                
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
            gui.myErrordlg(parentfig, ME.message, ME.identifier);
            return;
        end
        done = true;
    end

end
