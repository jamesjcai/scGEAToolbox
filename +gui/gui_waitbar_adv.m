function [f]=gui_waitbar_adv(f,p)
if nargin<2, p=1; end
if nargin<1 || isempty(f)
    f = waitbar(0,'Please wait...');
    pause(.5)
    fprintf('Processing your data...');
    waitbar(.67,f,'Processing your data');
    fprintf('... ');
    tic;
    return;
elseif isvalid(f) && strcmp(f.Tag,'TMWWaitbar')
    if p==1
        toc;    
        waitbar(1,f,'Finishing');
        pause(1);
        % fprintf('.......................done.\n');
        if isvalid(f)
            close(f);
        end
    else
        waitbar(p,f,'Processing your data');
    end
end
