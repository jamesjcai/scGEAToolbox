function [done] = callback_Harmonypy(src, ~)
done = false;

%[ok] = gui.i_confirmscript('Run Batch Integration (Harmony)?', ...
%    'py_harmonypy', 'python');
%if ~ok, return; end

[FigureHandle, sce, isui] = gui.gui_getfigsce(src);
if numel(unique(sce.c_batch_id)) < 2
    gui.myWarndlg(FigureHandle, 'No batch effect (all cells have the same SCE.C_BATCH_ID)');
    return;
end

extprogname = 'py_harmonypy';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
if isempty(wkdir), return; end
if ~gui.i_setpyenv, return; end



%{
usepylib=false;

answer = gui.myQuestdlg(FigureHandle, 'Using MATLAB engine for Python or Calling Python script?', ...
        'Engine Interface', ...
        'Use MATLAB Engine for Python','Call Python Script',...
            'Cancel','Use MATLAB Engine for Python');
            switch answer
                case 'Use MATLAB Engine for Python'
                        usepylib=true;
                    case 'Call Python Script'
                        usepylib=false;
                    case {'Cancel',''}
                        return;
                    otherwise
                        return;
                    end
                    %}

                    %fw=gui.gui_waitbar;

                    %try
                        
                        id = sce.c_batch_id;
                        if ~isnumeric(id)
                            id = grp2idx(sce.c_batch_id);
                            id = id(:);
                        end
                        [s] = run.py_harmonypy_new(sce.s, id, wkdir);
                        % [s] = run.py_harmonypy(sce.s, id);

                        if isempty(s) || isequal(sce.s, s)
                            % gui.gui_waitbar(fw);
                            errordlg("Harmonypy Running Error");
                            return;
                        end
                        sce.s = s;
                    % catch ME
                    %     %gui.gui_waitbar(fw,true);
                    %     errordlg(ME.message);
                    %     %rethrow(ME);
                    %     return;
                    % end
                    %gui.gui_waitbar(fw);

                    guidata(FigureHandle, sce);
                    done = true;
            end
