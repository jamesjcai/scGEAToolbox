function [done]=callback_Harmonypy(src,~)
    done=false;
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    if numel(unique(sce.c_batch_id))<2
        warndlg('No batch effect (SCE.C_BATCH_ID is empty)');
        return;
    end
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
        end
        fw=gui.gui_waitbar;
        try
            sce.s=run.harmonypy(sce.s,sce.c_batch_id);
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
        [file,path] = uigetfile('python.exe');
        if isequal(file,0)
           disp('User selected Cancel');
           return;
        else
           disp(['User selected ', fullfile(path,file)]);
           pyenv('Version',fullfile(path,file));
           done=true;
        end
end