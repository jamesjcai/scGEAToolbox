function [done, changed] = i_setwrkdir(preftagname, parentfig, askchange)
% I_SETWRKDIR - set working directory
%
% DONE is true when a working folder is set on return. CHANGED is true
% only when this call stored a new one.
%
% With ASKCHANGE false (the default), an already-set folder is accepted
% silently -- what the analyses that merely need a folder want. With
% ASKCHANGE true the current folder is shown and can be kept or replaced,
% which is what Setup > Set Working Folder is for: without it that menu
% could only ever set the folder once.
% see also: I_SETPYENV, I_SETRENV

if nargin < 3, askchange = false; end
if nargin < 2, parentfig = []; end
if nargin < 1, preftagname = 'externalwrkpath'; end
changed = false;
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end
[done] = false;

% ISPREF alone is not the question: the preference can exist and hold an
% empty value, which every caller then treats as "set up" and fails on.
% GUI.GUI_SETPRGMWKDIR errored outright in that state, so a stale empty
% preference took out every external tool that needs a working folder --
% Memento, CellBender, Monocle3, copykat, SCEVAN, DecontX, Seurat -- with
% no way to reach the dialog that would fix it. An empty value means not
% set up, and gets the same prompt as no value at all.
issetup = ispref('scgeatoolbox', preftagname) && ...
    ~isempty(getpref('scgeatoolbox', preftagname, []));

if ~issetup
    answer = gui.myQuestdlg(parentfig, ['Working directory has ' ...
        'not been set up. Locate a folder?']);
    if ~strcmp(answer, 'Yes'), return; end
    if ispc
        pathdefult = fullfile(getenv('USERPROFILE'), 'Documents');
    else
        pathdefult = '';
    end
    [done] = ix_setwdpath(pathdefult, parentfig);
    changed = done;
elseif askchange
    curpath = getpref('scgeatoolbox', preftagname);
    answer = gui.myQuestdlg(parentfig, curpath, 'Working Folder', ...
        {'Use this', 'Use another', 'Cancel'}, 'Use this');
    switch answer
        case 'Use this'
            done = true;
        case 'Use another'
            [done] = ix_setwdpath(curpath, parentfig);
            changed = done;
        otherwise
            % Cancel, or the dialog was closed: the folder stays set, but
            % the user did not confirm it, so report nothing done.
            return;
    end
else
    done = true;
end


function [done] = ix_setwdpath(deflt, parentfig)
        done = false;
        answer=gui.myQuestdlg(parentfig, 'Where to save working files?','',...
            {'Use Temporary Folder', ...
            'Select a Folder','Cancel'},'Use Temporary Folder');
        if isempty(answer), return; end
        switch answer
            case 'Select a Folder'
                [seltpath] = uigetdir(deflt);
                if pkg.i_isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure')
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
