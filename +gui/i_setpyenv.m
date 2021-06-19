function [done]=i_setpyenv
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