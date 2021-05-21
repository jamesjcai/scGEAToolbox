function [f]=gui_waitbar(f)
if nargin<1 || isempty(f)
    f = waitbar(0,'Please wait...');
    pause(.5)
    fprintf('Processing your data...');
    waitbar(.67,f,'Processing your data');    
    tic;
    return;
elseif isvalid(f) && strcmp(f.Tag,'TMWWaitbar')
    toc;    
    waitbar(1,f,'Finishing');
    pause(1);
    fprintf('...done.\n');
    if isvalid(f)
        close(f);
    end
end
