function [fw] = myWaitbar(parentfig, fw, witherror, mesg, newmesg, f)


if nargin < 6, f = []; end
if nargin < 5, newmesg = ''; end
if nargin < 4 || isempty(mesg), mesg = 'Processing your data'; end
if nargin < 3 || isempty(witherror), witherror = false; end
if nargin < 1, parentfig = []; end



    if ~gui.i_isuifig(parentfig)
        if nargin < 2 || isempty(fw)
            % hFig = gcf;
            hFig = get(groot,'CurrentFigure');
            fw = waitbar(0, 'Please wait...','Visible','off', ...
                'Units','pixels');
            try
                if ~isempty(hFig) && strcmp(get(hFig,'type'),'figure')
                    [~, newpos] = gui.i_getchildpos(hFig, fw);
                    fw.Position = newpos;
                end
            catch
            end
        
            fw.Visible = "on";
            pause(.5)
            fprintf('Processing your data...');
            fw = waitbar(0.618, fw, mesg);
            fprintf('... ');
            tic;
            return;
        elseif isvalid(fw) && strcmp(fw.Tag, 'TMWWaitbar') && ~isempty(newmesg) && isempty(f)
            newmesg = strrep(newmesg,'_','\_');
            fw = waitbar(.618, fw, newmesg);
        elseif isvalid(fw) && strcmp(fw.Tag, 'TMWWaitbar') && isempty(newmesg) && ~isempty(f)
            fw = waitbar(f, fw);
        elseif isvalid(fw) && strcmp(fw.Tag, 'TMWWaitbar') && ~isempty(newmesg) && ~isempty(f)
            newmesg = strrep(newmesg,'_','\_');
            fw = waitbar(f, fw, newmesg);
        elseif isvalid(fw) && strcmp(fw.Tag, 'TMWWaitbar')
            if ~witherror
                if nargin < 3 || isempty(mesg), mesg = 'Finishing'; end
                toc;
                fw = waitbar(1, fw, mesg);
                pause(1);
                % fprintf('.......................done.\n');
            end
            if isvalid(fw), close(fw); end
        end

    else

        if nargin < 2 || isempty(fw)
            % hFig = gcf;
            %hFig = get(groot,'CurrentFigure');
            %hUiFigHandle = findall(0, 'Type', 'figure', 'BeingDeleted', 'off');

            fw = uiprogressdlg(parentfig, 'Title', 'Please wait...', ...
                'Message', mesg);
            fprintf('Processing your data...');
            fw.Value = 0.618;
            fprintf('... ');
            tic;
            return;
        elseif isvalid(fw) && isa(fw, 'matlab.ui.dialog.ProgressDialog') && ...
                ~isempty(newmesg) && isempty(f)
            fw.Value = 0.618;
            fw.Message = newmesg;
        elseif isvalid(fw) && isa(fw, 'matlab.ui.dialog.ProgressDialog') && ...
                isempty(newmesg) && ~isempty(f)
            fw.Value = f;
        elseif isvalid(fw) && isa(fw, 'matlab.ui.dialog.ProgressDialog') && ...
                ~isempty(newmesg) && ~isempty(f)
            fw.Message = newmesg;
            fw.Value = f;
        elseif isvalid(fw) && isa(fw, 'matlab.ui.dialog.ProgressDialog')
            if ~witherror
                if nargin < 3 || isempty(mesg), mesg = 'Finishing'; end
                toc;
                fw.Value = 1;
                fw.Message = mesg;
                pause(1);
                % fprintf('.......................done.\n');
            end
            if isvalid(fw), close(fw); end
        end
    end
end

