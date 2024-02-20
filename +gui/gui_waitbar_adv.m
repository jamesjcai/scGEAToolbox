function [f] = gui_waitbar_adv(f, p, msg)
if nargin < 3, msg = []; end
if nargin < 2, p = 1; end
if nargin < 1 || isempty(f)

    hFig = gcf;
    %     if ~usejava('desktop')
    %         disp(msg);
    %         return;
    %     end
    f = waitbar(0, 'Please wait...','Visible','off','Units','pixels');

    [~, newpos] = gui.i_getchildpos(hFig, f);
    f.Position = newpos;
    f.Visible = "on";
    pause(.5)
    %     pause(.5)
    %     fprintf('Processing your data...');
    %     waitbar(.67,f,'Processing your data');
    %     fprintf('... ');
    tic;
    return;
elseif isvalid(f) && strcmp(f.Tag, 'TMWWaitbar')
    if p == 1
        toc;
        waitbar(1, f, 'Finishing');
        pause(1);
        % fprintf('.......................done.\n');
        if isvalid(f)
            close(f);
        end
    else
        if isempty(msg)
            msg = 'Processing your data';
        end
        %if ~usejava('desktop')
        %    disp(msg);
        %else
        msg = strrep(msg, '_', '\_');
        waitbar(p, f, msg);
        %end
    end
end
