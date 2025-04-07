function [done] = i_setpyenv(src, ~, parentfig)
% selpath = uigetdir;
%see also: I_SETRENV
[done] = false;

narginchk(2,3);
if nargin < 3
    [parentfig] = gui.gui_getfigsce(src);
end

x = pyenv;
if x.Version == "" %strlength(x.Executable)==0
    answer = gui.myQuestdlg(parentfig, ['Python environment ' ...
        'has not been set up. Locate python.exe?']);
    if strcmp('Yes', answer) 
        [done] = ix_setpyenv(x.Executable);
        if ~done
            return;
        end
        gui.myHelpdlg(parentfig, 'Python environment is set successfully.', '');
    else
        return;
    end
else
    answer = gui.myQuestdlg(parentfig, sprintf('%s', x.Executable), ...
        'Python Executable', ...
        {'Use this', 'Use another', 'Cancel'}, 'Use this');
    if isempty(answer), return; end
    switch answer
        case 'Use this'
            done = true;
        case 'Use another'
            if ~ix_setpyenv(x.Executable)
                return;
            end
            done = true;
            gui.myHelpdlg(parentfig, ...
            'Python environment is set successfully.');
        case {'Cancel', ''}
            return;
        otherwise
            return;
    end
end
if ~done
    gui.myWarndlg(parentfig, ...
        'Python environment is not set.', '');
end




    function [done] = ix_setpyenv(deflt)
    % selpath = uigetdir;
        done = false;
    
        if ispc
            [file, path] = uigetfile('python.exe', 'Select Python Interpreter', deflt);
        else
            [file, path] = uigetfile('python', 'Select Python Interpreter', deflt);
        end
        if isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure'), figure(parentfig); end
        if isequal(file, 0)
            return;
        else
            disp(['User selected: ', fullfile(path, file)]);
            try
                pyenv('Version', fullfile(path, file));
            catch ME
                content = regexprep(ME.message, '<.*?>', '' ) ;
                gui.myErrordlg(parentfig, content);
                return;
            end
            done = true;
        end
    end

end
