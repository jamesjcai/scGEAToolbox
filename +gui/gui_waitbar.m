function [f]=gui_waitbar(f,witherror,mesg)

if nargin<3 || isempty(mesg), mesg='Processing your data'; end
if nargin<2 || isempty(witherror), witherror=false; end
if nargin<1 || isempty(f)    
    if ~usejava('desktop')
        disp(mesg);
        f=[];
    else
        f = waitbar(0,'Please wait...');
        pause(.5)
        fprintf('Processing your data...');
        waitbar(.67,f,mesg);
        fprintf('... ');
    end
    tic;    
    return;
elseif isvalid(f) && strcmp(f.Tag,'TMWWaitbar')
    if ~witherror
        if nargin<3 || isempty(mesg), mesg='Finishing'; end
        toc;
        waitbar(1,f,mesg);
        pause(1);
        % fprintf('.......................done.\n');
    end
    if isvalid(f), close(f); end
end
