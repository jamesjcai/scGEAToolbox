function [isDoublet, doubletscore, methodtag, done] = callback_DoubletDetection(src, ~)

done = false;
isDoublet = [];
doubletscore = [];
methodtag = [];

%[ok] = gui.i_confirmscript('Run Detect Doublets (Scrublet)?', ...
%    'py_scrublet', 'python');
%if ~ok, return; end

[FigureHandle, sce] = gui.gui_getfigsce(src);


extprogname = 'py_scrublet';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
if isempty(wkdir), return; end
if ~gui.i_setpyenv, return; end


% methodtag = 'scrublet';

% methodtag=questdlg('Which method?','',...
%     'scrublet','doubletdetection','scrublet');

%fw=gui.gui_waitbar;
% try
%     switch methodtag
%         case 'scrublet'
try
    
[isDoublet, doubletscore] = run.py_scrublet_new(sce.X, wkdir);
            %                 case 'doubletdetection'
            %                     [isDoublet,doubletscore]=run.py_doubletdetection(sce.X);
    %     otherwise
    %         return;
    % end
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
