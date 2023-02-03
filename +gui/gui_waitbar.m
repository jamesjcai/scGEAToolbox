function [f]=gui_waitbar(f,witherror,mesg)

if nargin<3, mesg='Processing your data'; end
if nargin<2 || isempty(witherror), witherror=false; end
if nargin<1 || isempty(f)
    f = waitbar(0,'Please wait...');
    pause(.5)
    fprintf('Processing your data...');
    waitbar(.67,f,mesg);
    fprintf('... ');
    tic;
    return;
elseif isvalid(f) && strcmp(f.Tag,'TMWWaitbar')
    if ~witherror
        toc;
        waitbar(1,f,'Finishing');
        pause(1);
        % fprintf('.......................done.\n');
    end
    if isvalid(f), close(f); end
end
