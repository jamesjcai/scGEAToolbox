function [done] = i_setrenv(src, ~)

%see also: I_SETPYENV, I_SETEXTWD
[done] = false;
if nargin<1
    parentfig = [];
else
    [parentfig, ~] = gui.gui_getfigsce(src);
end
preftagname = 'rexecutablepath';

if ismac
    if ispref('scgeatoolbox', preftagname)
        s = getpref('scgeatoolbox', preftagname);
    else
        s = '/Library/Frameworks/R.framework/Versions/4.4/Resources/bin/R';
    end
    
if gui.i_isuifig(parentfig)
    new_s = gui.myInputdlg({'Path to R Executable (e.g., /Library/Frameworks/R.framework/Versions/4.4/Resources/bin/R)'}, ...
                           'Set up R Environment', {s}, parentfig);
else
    new_s = inputdlg('Path to R Executable (e.g., /Library/Frameworks/R.framework/Versions/4.4/Resources/bin/R)', ...
                     'Set up R Environment', [1 100], {s});
end    


    if isempty(new_s), return; end
    if ~isempty(new_s{1}) && strlength(new_s{1})
        setpref('scgeatoolbox', preftagname, new_s{1});
        done = true;
    end
    return;
end

if ~ispref('scgeatoolbox', preftagname)
    answer = gui.myQuestdlg(parentfig, 'R environment has not been set up. Locate R executable Rscript.exe?');
    if ~strcmp(answer, 'Yes'), return; end
    if ispc
        rpathdefult = pkg.i_findrpath;
        if iscell(rpathdefult)
            rpathdefult = rpathdefult{1};
        end
    else
        rpathdefult = '';
    end
    [done] = ix_setrenv(rpathdefult);
    % setpref('scgeatoolbox',preftagname,s);
else
    s = getpref('scgeatoolbox', preftagname);
    if isempty(s)
        rmpref('scgeatoolbox', preftagname);
        done = false;
    else
        answer = gui.myQuestdlg(parentfig, sprintf('%s', s), ...
            'Path to R Executable', ...
            {'Use this', 'Use another', 'Cancel'}, 'Use this');
        switch answer
            case 'Use this'
                done = true;
            case 'Use another'
                [done] = ix_setrenv(s);
        end
    end
end
if done
    Rpath = getpref('scgeatoolbox', preftagname);
    if ispc
        Rexec = fullfile(Rpath, 'Rscript.exe');
    else
        Rexec = fullfile(Rpath, 'Rscript');
    end
    if exist(Rexec, 'file')
        gui.myHelpdlg(parentfig, "R environment is set successfully.");
    else
        if ispc
            gui.myErrordlg(parentfig, "R environment is set with error.", '');
            done = false;
        else
            Rexec = fullfile(Rpath, 'R');
            if exist(Rexec, 'file')
                gui.myHelpdlg(parentfig, "R environment is set successfully.", '');
            else
                gui.myErrordlg(parentfig, "R environment is set with error.", '');
                done = false;
            end
        end
    end
end

    function [done] = ix_setrenv(deflt)
        % selpath = uigetdir;
        done = false;
        if ispc
            [file, path] = uigetfile('Rscript.exe', 'Select R Interpreter', deflt);
        else
            [file, path] = uigetfile('Rscript', 'Select R Interpreter', deflt);
        end
        if isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure'), figure(parentfig); end
        if isequal(file, 0) 
            return; 
        else
            disp(['User selected: ', fullfile(path, file)]);
            % fullfile(path)
            try
                setpref('scgeatoolbox', preftagname, fullfile(path));
            catch ME
                gui.myErrordlg(parentfig, ME.message, ME.identifier);
                return;
            end
            done = true;
        end
    end

end
