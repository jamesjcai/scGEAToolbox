function [isDoublet, doubletscore, methodtag, done] = callback_DoubletDetection(src, ~)

done = false;
isDoublet = [];
doubletscore = [];
methodtag = [];

[ok] = gui.i_confirmscript('Run Detect Doublets (Scrublet)?', ...
    'py_scrublet', 'python');
if ~ok, return; end

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

if ~gui.i_setpyenv, return; end

methodtag = 'scrublet';

% methodtag=questdlg('Which method?','',...
%     'scrublet','doubletdetection','scrublet');

%fw=gui.gui_waitbar;
try
    switch methodtag
        case 'scrublet'
            [isDoublet, doubletscore] = run.py_scrublet(sce.X);
            %                 case 'doubletdetection'
            %                     [isDoublet,doubletscore]=run.py_doubletdetection(sce.X);
        otherwise
            return;
    end
    if isempty(isDoublet) || isempty(doubletscore)
        %gui.gui_waitbar(fw);
        errordlg("Running Error.");
        return;
    end
catch ME
    %gui.gui_waitbar(fw);
    errordlg(ME.message);
    rethrow(ME);
end
%gui.gui_waitbar(fw);
guidata(FigureHandle, sce);
done = true;
end
