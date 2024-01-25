function [f] = uiwaitbar(fig, f, witherror, mesg, newmesg)
if nargin < 5, newmesg = ''; end
if nargin < 4 || isempty(mesg), mesg = 'Processing your data'; end
if nargin < 3 || isempty(witherror), witherror = false; end
if nargin < 2 || isempty(f)
    f = uiprogressdlg(fig, 'title','Please wait...',...
        'Message','Processing the data');
    pause(.5)
    fprintf('Processing your data...');
    f.Value = .618;
    fprintf('... ');
    tic;
    return;
elseif isvalid(f) && ~isempty(newmesg)
    f.Message = newmesg;
elseif isvalid(f)
    if ~witherror
        if nargin < 3 || isempty(mesg), mesg = 'Finishing'; end
        toc;
        f.Value =1;
        f.Message = mesg;
        pause(1);        
    end
    if isvalid(f), close(f); end
end
