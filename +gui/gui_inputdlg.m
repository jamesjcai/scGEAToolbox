function [numfig]=gui_inputdlg
    numfig=10;
    prompt = {'Enter number of genes (1-50):'};
    dlgtitle = 'Number of marker genes';
    % answer = inputdlg(prompt,dlgtitle,[1 50],{'10'});
    
[answer] = timeoutDlg(@inputdlg,10,...
            prompt,dlgtitle,[1 50],{'10'});
if iscell(answer)
    answer=answer{1};
end
    if isempty(answer)
        return;
    else
        try
            numfig=str2double(answer);
        catch ME
            errordlg(ME.message);
            return;
        end
    end
    if ~(numfig>0 && numfig<=50)
        errordlg('Invalid number of figures');
        return;
    end
end


function varargout = timeoutDlg(dlg, delay, varargin)
% Dialog function with timeout property
% dlg is a handle to the dialog function to be called
% delay is the length of the delay in seconds
% other input arguments as required by the dialog
% EXAMPLE FUNCTION-CALL
% To display an input dialog box (REFER MATLAB HELP DOC) with a 
% timeout = 6 second say, the function call would be:
%
% [matrix_size_value, colormap_string] = timeoutdlg(@inputdlg, 6, ...
%                                {'Enter matrix size:','Enter colormap name:'}, ...
%                                'Input for peaks function', 1, {'20','hsv'})
% Setup a timer to close the dialog in a moment
f1 = findall(0, 'Type', 'figures');
t = timer('TimerFcn', {@closeit f1}, 'StartDelay', delay);
start(t);
% Call the dialog
retvals = dlg(varargin{:});
if numel(retvals) == nargout
      varargout = retvals(:);
else
      varargout = cell(1, nargout);
end
% Delete the timer
if strcmp(t.Running, 'on')
       stop(t);
end
delete(t);
    function closeit(src, event, f1)
    %disp('Time out!');
    disp('Using defult value = 30');
    f2 = findall(0, 'Type', 'figure');
    fnew = setdiff(f2, f1);
    if ishandle(fnew)
          close(fnew);
    end
    end
end
