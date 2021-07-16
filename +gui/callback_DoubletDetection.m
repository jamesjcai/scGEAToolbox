function [isDoublet,doubletscore,methodtag,done]=callback_DoubletDetection(src,~)
    done=false;
    isDoublet=[];
    doubletscore=[];
    methodtag=[];
    
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

if ~gui.i_setpyenv, return; end

methodtag='scrublet';

% methodtag=questdlg('Which method?','',...
%     'scrublet','doubletdetection','scrublet');
        
        fw=gui.gui_waitbar;
        try
            switch methodtag
                case 'scrublet'
                    [isDoublet,doubletscore]=run.scrublet(sce.X);
%                 case 'doubletdetection'
%                     [isDoublet,doubletscore]=run.doubletdetection(sce.X);
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
   guidata(FigureHandle,sce);
   done=true;
end
