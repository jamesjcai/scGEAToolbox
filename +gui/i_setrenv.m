function [done] = i_setrenv(~, ~)

%see also: I_SETPYENV
[done] = false;

% Rpath=pkg.FindRpath;
% if isempty(Rpath)
%     warndlg('R is not installed.','');
% else
%     if iscell(Rpath)
%         s=Rpath{1};
%         for k=2:length(Rpath)
%             s=sprintf('%s\n%s',s,Rpath{k});
%         end
%         s=sprintf('%s %s',s,' (default)');
%     else
%         s=Rpath;
%     end

%
if ~ispref('scgeatoolbox', 'rexecutablepath')
    answer = questdlg('Select Path to R Executable', '');
    if ~strcmp(answer, 'Yes'), return; end
    if ispc
        rpathdefult = pkg.FindRpath;
        if iscell(rpathdefult)
            rpathdefult = rpathdefult{1};
        end
    else
        rpathdefult = '';
    end
    [done] = ix_setrenv(rpathdefult);
    % setpref('scgeatoolbox','rexecutablepath',s);
else
    s = getpref('scgeatoolbox', 'rexecutablepath');
    answer = questdlg(sprintf('%s', s), ...
        'Path to R Executable', ...
        'Use this', 'Use another', 'Cancel', 'Use this');
    switch answer
        case 'Use this'
            done = true;
        case 'Use another'
            [done] = ix_setrenv(s);
    end
end
if done
    Rpath = getpref('scgeatoolbox', 'rexecutablepath');
    if ispc
        Rexec = fullfile(Rpath, 'Rscript.exe');
    else
        Rexec = fullfile(Rpath, 'Rscript');
    end
    if exist(Rexec, 'file')
        helpdlg("R environment is set successfully.", '');
    else
        if ispc
            errordlg("R environment is set with error.", '');
            done = false;
        else
            Rexec = fullfile(Rpath, 'R');
            if exist(Rexec, 'file')
                helpdlg("R environment is set successfully.", '');
            else
                errordlg("R environment is set with error.", '');
                done = false;
            end
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
if isequal(file, 0)
    %disp('User selected Cancel');
    return;
else
    disp(['User selected: ', fullfile(path, file)]);
    % fullfile(path)
    try
        setpref('scgeatoolbox', 'rexecutablepath', fullfile(path));
    catch ME
        errordlg(ME.message);
        return;
    end
    done = true;
end
end