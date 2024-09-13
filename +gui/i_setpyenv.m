function [done] = i_setpyenv(~, ~)
% selpath = uigetdir;
%see also: I_SETRENV
[done] = false;

x = pyenv;
if x.Version == "" %strlength(x.Executable)==0
    answer = questdlg('Python environment has not been set up. Locate python.exe?');
    if strcmp('Yes', answer) 
        [done] = ix_setpyenv(x.Executable);
        if ~done
            return;
        end
        waitfor(helpdlg('Python environment is set successfully.', ''));
    else
        return;
    end
else
    answer = questdlg(sprintf('%s', x.Executable), ...
        'Python Executable', ...
        'Use this', 'Use another', 'Cancel', 'Use this');
    switch answer
        case 'Use this'
            done = true;
        case 'Use another'
            if ~ix_setpyenv(x.Executable)
                return;
            end
            done = true;
            waitfor(helpdlg('Python environment is set successfully.', ''));
        case {'Cancel', ''}
            return;
        otherwise
            return;
    end
end
if ~done
    warndlg('Python environment is not set.', '');
end

end


    function [done] = ix_setpyenv(deflt)
    % selpath = uigetdir;
        done = false;
    
        if ispc
            [file, path] = uigetfile('python.exe', 'Select Python Interpreter', deflt);
        else
            [file, path] = uigetfile('python', 'Select Python Interpreter', deflt);
        end
        if isequal(file, 0)
            %disp('User selected Cancel');
            return;
        else
            disp(['User selected: ', fullfile(path, file)]);
            try
                pyenv('Version', fullfile(path, file));
            catch ME
                errordlg(ME.message);
                return;
            end
            done = true;
        end
    end
