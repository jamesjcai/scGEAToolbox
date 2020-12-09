function [f]=gui_waitbar(f)
if nargin<1 || isempty(f)    
    f = waitbar(0,'Please wait...');
    pause(.5)
    waitbar(.67,f,'Processing your data');
    return;
elseif strcmp(f.Tag,'TMWWaitbar')
    waitbar(1,f,'Finishing');
    pause(1); 
    close(f);    
end
