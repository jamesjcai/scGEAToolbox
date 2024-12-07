function [done] = i_setscimilaritymodelpath(~, ~)
% selpath = uigetdir;
%see also: I_SETRENV
[done] = false;

preftagname = 'scimilmodelpath';
if ispref('scgeatoolbox', preftagname)
    x = getpref('scgeatoolbox', preftagname);
else
    x = '';
end

if isempty(x) || ~isfolder(x)
    answer = questdlg('Scimilarity model path has not been set up. Locate it?');
    if strcmp('Yes', answer) 
        [done] = ix_setpyenv;
        if ~done, return; end
        waitfor(helpdlg('Scimilarity model path is set successfully.', ''));
    else
        return;
    end
else
    answer = questdlg(sprintf('%s', x), ...
        'Model Path', ...
        'Use this', 'Use another', 'Cancel', 'Use this');
    switch answer
        case 'Use this'
            done = true;
        case 'Use another'
            if ~ix_setpyenv
                return;
            end
            done = true;
            waitfor(helpdlg('Scimilarity model path is set successfully.', ''));
        case {'Cancel', ''}
            return;
        otherwise
            return;
    end
end
if ~done
    warndlg('Scimilarity model path is not set.', '');
end

end


    function [done] = ix_setpyenv(deflt)
    % selpath = uigetdir;
        done = false;
    
        if ispc
            [file, path] = uigetfile('python.exe', 'Select Python Interpreter', deflt);
        else
            [file, path] = uigetfile('python', 'Select Python Interpreter', deflt);
        end
        if isequal(file, 0)
            %disp('User selected Cancel');
            return;
        else
            disp(['User selected: ', fullfile(path, file)]);
            try
                pyenv('Version', fullfile(path, file));
            catch ME
                content = regexprep(ME.message, '<.*?>', '' ) ;
                waitfor(errordlg(content,''));
                return;
            end
            done = true;
        end
    end
