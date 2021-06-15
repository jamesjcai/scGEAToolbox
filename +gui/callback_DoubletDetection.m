function [isDoublet,doubletscore,done]=callback_DoubletDetection(src,~)
    done=false;
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    x=pyenv;
    if isempty(x.Executable)
        i_setpyenv
    else
        answer = questdlg(sprintf('%s',x.Executable), ...
            'Python Executable', ...
            'Use this','Use another','Cancel','Use this');        
        switch answer
            case 'Use this'
            case 'Use another'
                if ~i_setpyenv, return; end                    
            case {'Cancel',''}
                return;
            otherwise
                return;
        end

        
        fw=gui.gui_waitbar;
        try
            [isDoublet,doubletscore]=run.doubletdetection(sce.X);
            if isempty(isDoublet)
                gui.gui_waitbar(fw);
                errordlg("doubletdetection Running Error");
                return;
            end            
        catch ME
            gui.gui_waitbar(fw);
            errordlg(ME.message);
            rethrow(ME);
        end 
            gui.gui_waitbar(fw);        
    end
   guidata(FigureHandle,sce);
   done=true;
end


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