function [isDoublet,doubletscore,methodtag,done]=callback_DoubletDetection(src,~)
    done=false;
    isDoublet=[];
    doubletscore=[];
    methodtag=[];
    
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
                if ~gui.i_setpyenv, return; end                    
            case {'Cancel',''}
                return;
            otherwise
                return;
        end

methodtag=questdlg('Which method?','',...
    'scrublet','doubletdetection','scrublet');
        
        fw=gui.gui_waitbar;
        try
            switch methodtag
                case 'scrublet'
                    [isDoublet,doubletscore]=run.scrublet(sce.X);
                case 'doubletdetection'
                    [isDoublet,doubletscore]=run.doubletdetection(sce.X);
                otherwise
                    return;
            end
            if isempty(isDoublet)||isempty(doubletscore)
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
