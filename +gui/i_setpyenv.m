function [done]=i_setpyenv
        % selpath = uigetdir;
           
    x=pyenv;
    if isempty(x.Executable)
        [done]=ix_setpyenv;
    else
        [done]=false;
        answer = questdlg(sprintf('%s',x.Executable), ...
            'Python Executable', ...
            'Use this','Use another','Cancel','Use this');        
        switch answer
            case 'Use this'
            case 'Use another'
                if ~ix_setpyenv, return; end                    
            case {'Cancel',''}
                return;
            otherwise
                return;
        end
        [done]=true;
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
           pyenv('Version',fullfile(path,file));
           done=true;
        end
end