function [f] = gui_waitbar_simple(f, witherror, mesg, newmesg)
if nargin < 4, newmesg = ''; end
if nargin < 3 || isempty(mesg), mesg = 'Processing your data'; end
if nargin < 2 || isempty(witherror), witherror = false; end
if nargin < 1 || isempty(f)
    
    %hFig = gcf;
    f = waitbar(0, 'Please wait...','Visible','off','Units','pixels');
    %[~, newpos] = gui.i_getchildpos(hFig, f);
    %f.Position = newpos;

    f.Visible = "on";
    pause(.5)
    fprintf('Processing your data...');
    waitbar(.618, f, mesg);
    fprintf('... ');
    tic;
    return;
elseif isvalid(f) && strcmp(f.Tag, 'TMWWaitbar') && ~isempty(newmesg)
    waitbar(.618, f, newmesg);
elseif isvalid(f) && strcmp(f.Tag, 'TMWWaitbar')
    if ~witherror
        if nargin < 3 || isempty(mesg), mesg = 'Finishing'; end
        toc;
        waitbar(1, f, mesg);
        pause(1);
        % fprintf('.......................done.\n');
    end
    if isvalid(f), close(f); end
end
