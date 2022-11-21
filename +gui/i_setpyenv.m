function [done]=i_setpyenv(~,~)
        % selpath = uigetdir;
%see also: I_SETRENV
        [done]=false;

    x=pyenv;
    if strlength(x.Executable)==0
        answer=questdlg('Python environment has not been set up. Locate python.exe?');
        if strcmp('Yes',answer)
            [done]=ix_setpyenv;
        else
            return;
        end
    else
        answer = questdlg(sprintf('%s',x.Executable), ...
            'Python Executable', ...
            'Use this','Use another','Cancel','Use this');        
        switch answer
            case 'Use this'
                done=true;
            case 'Use another'
                if ~ix_setpyenv
                    return;
                else
                    done=true;
                end
            case {'Cancel',''}
                return;
            otherwise
                return;
        end
    end
    if ~done
        warndlg('Python environment is not set.','');
    end
    
end


function [done]=ix_setpyenv
        % selpath = uigetdir;
        done=false;        
        
        if ispc
            [file,path] = uigetfile('python.exe','Select Python Interpreter');
        else
            [file,path] = uigetfile('python','Select Python Interpreter');
        end
        if isequal(file,0)
           %disp('User selected Cancel');
           return;
        else
           disp(['User selected: ', fullfile(path,file)]);
           try
                pyenv('Version',fullfile(path,file));
           catch ME
               errordlg(ME.message);
               return;
           end
           done=true;
        end
end