classdef myFigure < handle
    % MYFIGURE A wrapper class to create customized MATLAB figures with toolbar
    %
    % Usage:
    %   hx = gui.myFigure(parentfig)           % auto-detect: uifigure
    %                                           % parent -> uifigure child
    %   hx = gui.myFigure(parentfig, true)      % force classic figure
    %
    % When parentfig is a uifigure, the constructor creates a uifigure
    % with uiaxes. This is suitable for simple single-axes plots.
    %
    % However, classic figure functions such as subplot, datacursormode,
    % ginput, gca, and hold do NOT work on uifigures. If your code uses
    % any of these, pass forceclassic=true to create a classic figure
    % instead. Theme from the parent uifigure is still inherited.
    %
    % Examples:
    %   % Simple single-axes plot (uifigure OK):
    %   hx = gui.myFigure(parentfig);
    %   plot(hx.AxHandle, x, y);
    %   hx.show(parentfig);
    %
    %   % Multi-panel plot needing subplot (must force classic):
    %   hx = gui.myFigure(parentfig, true);
    %   hAx1 = subplot(2,1,1);
    %   hAx2 = subplot(2,1,2);
    %   hx.show(parentfig);
    %
    % Drawing into an existing window:
    %   A button on a myFigure window that shows another view should not
    %   open a third window over the main app and this one. Wrap the call
    %   that builds the view in DRAWINTO instead:
    %
    %   gui.myFigure.drawInto(hFig, @() gui.i_dotplot(X, g, c, cL, genes));
    %
    %   The first myFigure that the function creates takes over hFig rather
    %   than opening a window. The view that was there is kept, hidden, and a
    %   Back button on the new view brings it back as it was left. Pass
    %   KeepCurrent=false to discard the current view instead, for a pair
    %   of views that already switch to each other.

    properties
        FigHandle % Handle to the MATLAB figure
        AxHandle
    end

    properties (SetAccess = private)
        % The standard "Enlarge Window" toolbar button (GUI.I_TOGGLEWINSIZE).
        % Exposed so a plot whose own layout depends on the window size can
        % replace its ClickedCallback with one that wraps the toggle.
        SizeButton
    end

    properties (Access = private)
        tb
        tb2
        tbv = cell(13, 1) % Toolbar button handles
        % True when this view was drawn into an existing window by DRAWINTO.
        % SHOW then leaves the window where the user has it.
        Adopted = false
        % Standard buttons hidden while another view is drawn over this one.
        SuspendedButtons = gobjects(0)
    end

    properties (Constant, Access = private)
        % Figure properties a view may set for itself. They are saved with
        % a view that is put aside and reset for the view drawn over it.
        ViewProps = {'Name', 'Colormap', 'Pointer', 'WindowButtonDownFcn', ...
            'WindowButtonMotionFcn', 'WindowButtonUpFcn', 'WindowScrollWheelFcn', ...
            'WindowKeyPressFcn', 'WindowKeyReleaseFcn', 'KeyPressFcn', ...
            'KeyReleaseFcn', 'SizeChangedFcn', 'ThemeChangedFcn'}
    end

    methods
        % Constructor
        function obj = myFigure(parentfig, forceclassic)
            if nargin<2, forceclassic = true; end
            if nargin<1, parentfig = []; end

            useuifig = ~forceclassic && gui.i_isuifig(parentfig);
            pending = gui.myFigure.i_takehost(useuifig);
            if ~isempty(pending)
                obj.i_adopt(pending);
            elseif useuifig
                obj.FigHandle = uifigure('Name', '', ...
                    'Visible',"off","ToolBar","none");

                if ~isMATLABReleaseOlderThan('R2025a')
                    try
                        theme(obj.FigHandle, parentfig.Theme.BaseColorStyle);
                    catch
                        % theme() may not exist or parent has no Theme; skip styling
                    end
                end
            else
                obj.FigHandle = figure('Name', '', ...
                    'NumberTitle', 'on', 'Visible',"off", ...
                    'ToolBar','figure', ...
                    'MenuBar','figure', ...
                    'WindowStyle','normal', ...
                    'Position', [1, 1, 560, 420]);

                % --- Apply theme if available ---
                if ~isempty(parentfig) && isprop(parentfig,'Theme') ...
                        && ~isMATLABReleaseOlderThan('R2025a')
                    try
                        theme(obj.FigHandle, parentfig.Theme.BaseColorStyle);
                    catch
                        % Ignore if theme fails
                    end
                end
            end

            if useuifig
                obj.tb = uitoolbar(obj.FigHandle);
                obj.AxHandle = uiaxes(obj.FigHandle);
                obj.AxHandle.Position = [50, 30, obj.FigHandle.Position(3)-100, ...
                    obj.FigHandle.Position(4)-60];
                obj.i_addstandardbuttons('off');
            else
                obj.AxHandle = axes('Parent', obj.FigHandle);
                % The figure's own toolbar, by tag: a window this view was
                % drawn into may still hold the hidden toolbar of the view
                % underneath it.
                obj.tb = findall(obj.FigHandle, 'Tag', 'FigureToolBar');
                if isempty(obj.tb)
                    obj.tb = uitoolbar(obj.FigHandle);
                end
                % Remove undesirable default tools if present
                delete(findall(obj.FigHandle, 'Tag', 'DataManager.Linking'));
                delete(findall(obj.FigHandle, 'Tag', 'Standard.OpenInspector'));
                obj.i_addstandardbuttons('on');
            end

            if obj.Adopted && ~isempty(getappdata(obj.FigHandle, 'myFigureStack'))
                % First on the custom toolbar, ahead of the view's own buttons.
                host = obj.FigHandle;
                obj.addCustomButton('off', @(~, ~) gui.myFigure.goBack(host), ...
                    gui.myFigure.backIcon(), 'Back to previous view');
            end
            setappdata(obj.FigHandle, 'myFigureObj', obj);
        end

        function in_darkmode(obj, ~, ~)
            if isprop(obj.FigHandle,'Theme') ...
                    && ~isMATLABReleaseOlderThan('R2025a')
                try
                    if strcmp('light', obj.FigHandle.Theme.BaseColorStyle)
                        theme(obj.FigHandle, 'dark');
                    else
                        theme(obj.FigHandle, 'light');
                    end
                catch
                    % Ignore if theme fails
                end
            end
        end

        function ptvshow(obj, flag)
            obj.i_setprop(flag, 'Visible', 'on', 'off');
        end

        function ptvenable(obj, flag)
            obj.i_setprop(flag, 'Enable', 'on', 'off');
        end

        function centerto(obj, parentfig)
            gui.i_movegui2parent(obj.FigHandle, parentfig);
        end

        function show(obj, parentfig)
            if nargin<2, parentfig = []; end
            if obj.Adopted
                % The window is already where the user put it.
                obj.FigHandle.Visible = true;
                figure(obj.FigHandle);
                return
            end
            centerto(obj, parentfig);
            obj.FigHandle.Visible = true;
        end

        function pt = addCustomButton(obj, sepTag, callback, imgFil, tooltipTxt)
            % Returns the button so a caller can change its icon or tooltip.
            if isempty(obj.tb2)
                obj.tb2 = uitoolbar(obj.FigHandle);
            end
            pt = gui.i_addbutton2fig(obj.tb2, sepTag, callback, imgFil, tooltipTxt);
        end

        function setTitle(obj, titleStr)
            obj.FigHandle.Name = titleStr;
        end

        function closeFigure(obj)
            if pkg.i_isvalid(obj.FigHandle)
                close(obj.FigHandle);
            end
        end

        %% Destructor
        function delete(obj)
            closeFigure(obj);
        end
    end

    methods (Static)
        function drawInto(hostfig, drawfcn, options)
            % DRAWINTO Draw a view into an existing myFigure window.
            %   gui.myFigure.drawInto(hostfig, drawfcn) calls DRAWFCN, and
            %   the first gui.myFigure it creates is drawn into HOSTFIG
            %   instead of a new window. The view already in HOSTFIG is put
            %   aside, and the new view gets a Back button that restores it.
            %
            %   gui.myFigure.drawInto(..., KeepCurrent=false) discards the
            %   current view; no Back button leads to it.
            %
            %   HOSTFIG must be a window made by gui.myFigure, of the same
            %   kind (figure or uifigure) as the one DRAWFCN asks for.
            %   Otherwise DRAWFCN opens its own window, as it would unwrapped.
            arguments
                hostfig
                drawfcn (1,1) function_handle
                options.KeepCurrent (1,1) logical = true
            end
            if pkg.i_isvalid(hostfig) && isappdata(hostfig, 'myFigureObj')
                setappdata(groot, 'myFigureHost', ...
                    struct('Fig', hostfig, 'Keep', options.KeepCurrent));
                % Unclaimed when DRAWFCN is cancelled or fails before it
                % creates a figure; the next unrelated myFigure must not
                % land in this window.
                cleanupObj = onCleanup(@() gui.myFigure.i_clearhost());
            end
            drawfcn();
        end

        function cdata = backIcon()
            % BACKICON The Back button's icon, a mirrored green arrow.
            persistent icon
            if isempty(icon)
                imgPath = fullfile(fileparts(mfilename('fullpath')), '..', ...
                    'assets', 'Images', 'greenarrowicon.gif');
                try
                    [im, map] = imread(imgPath);
                    icon = fliplr(ind2rgb(im, map));
                catch
                    % icon file missing; a plain button still works
                    icon = zeros(16, 16, 3);
                end
            end
            cdata = icon;
        end

        function goBack(hostfig)
            % GOBACK Replace the current view with the one it was drawn over.
            stack = getappdata(hostfig, 'myFigureStack');
            if isempty(stack), return; end
            saved = stack(end);
            stack(end) = [];

            current = getappdata(hostfig, 'myFigureObj');
            if ~isempty(current) && isvalid(current)
                current.i_release();
            end
            delete(hostfig.Children);

            kids = saved.Panel.Children;
            for k = numel(kids):-1:1
                kids(k).Parent = hostfig;
            end
            delete(saved.Panel);
            hidden = saved.Hidden(isvalid(saved.Hidden));
            set(hidden, 'Visible', 'on', 'HandleVisibility', 'on');
            set(hostfig, saved.FigProps);
            if pkg.i_isvalid(saved.CurrentAxes)
                hostfig.CurrentAxes = saved.CurrentAxes;
            end

            saved.Obj.i_resume();
            setappdata(hostfig, 'myFigureObj', saved.Obj);
            setappdata(hostfig, 'myFigureStack', stack);
            figure(hostfig);
        end
    end

    methods (Access = private)
        function i_setprop(obj, flag, prop, onVal, offVal)
            von = obj.tbv(flag);
            von = von(~cellfun(@isempty, von));
            voff = obj.tbv(~flag);
            voff = voff(~cellfun(@isempty, voff));
            if ~isempty(von), set([von{:}], prop, onVal); end
            if ~isempty(voff), set([voff{:}], prop, offVal); end
        end

        function i_addstandardbuttons(obj, firstSep)
            obj.tbv{1} = gui.i_addbutton2fig(obj.tb, firstSep, @gui.i_invertcolor, 'INVERT.gif', 'Invert Colors');
            obj.tbv{2} = gui.i_addbutton2fig(obj.tb, 'off', @gui.i_linksubplots, "keyframes-minus.jpg", "Link Subplots");
            obj.tbv{3} = gui.i_addbutton2fig(obj.tb, 'off', {@gui.i_setboxon, obj.FigHandle}, 'border-out.jpg', 'Box ON/OFF');
            obj.tbv{4} = gui.i_addbutton2fig(obj.tb, 'off', @gui.i_renametitle, "align-top-box.jpg", 'Add/Edit Title');
            obj.tbv{5} = gui.i_addbutton2fig(obj.tb, 'off', {@gui.i_pickcolor, false}, 'color-wheel.jpg', 'Pick a New Colormap...');
            obj.tbv{6} = gui.i_addbutton2fig(obj.tb, 'off', @gui.i_changefontsize, 'text-size.jpg', 'Change Font Size');
            % The two size buttons sit right after Change Font Size: a toolbar
            % clips from the right when the window is narrow, which is exactly
            % when they are needed. Less used tools (export, camera, theme) follow.
            obj.tbv{11} = gui.i_addbutton2fig(obj.tb, 'on', {@gui.i_resizewin, obj.FigHandle}, 'scale-frame-reduce.jpg', 'Resize Plot Window');
            obj.tbv{13} = gui.i_addbutton2fig(obj.tb, 'off', {@gui.i_togglewinsize, obj.FigHandle}, 'scale-frame-enlarge.jpg', 'Enlarge Window');
            obj.SizeButton = obj.tbv{13};
            obj.tbv{7} = gui.i_addbutton2fig(obj.tb, 'on', {@gui.i_savemainfig, 3}, "presentation.jpg", 'Save Figure to PowerPoint File...');
            obj.tbv{8} = gui.i_addbutton2fig(obj.tb, 'off', {@gui.i_savemainfig, 2, obj.FigHandle, obj.AxHandle}, "jpg-format.jpg", 'Save Figure as Graphic File...');
            obj.tbv{9} = gui.i_addbutton2fig(obj.tb, 'off', {@gui.i_savemainfig, 1, obj.FigHandle, obj.AxHandle}, "svg-format.jpg", 'Save Figure as SVG File...');
            obj.tbv{10} = gui.i_3dcamera(obj.tb, '', false, obj.FigHandle, obj.AxHandle);
            obj.tbv{12} = gui.i_addbutton2fig(obj.tb, 'on', @obj.in_darkmode, 'demoIcon.gif', 'Light/Dark Mode');
        end

        function i_adopt(obj, pending)
            host = pending.Fig;
            previous = getappdata(host, 'myFigureObj');
            stack = getappdata(host, 'myFigureStack');
            if isempty(stack)
                stack = struct('Obj', {}, 'Panel', {}, 'Hidden', {}, ...
                    'FigProps', {}, 'CurrentAxes', {});
            end

            if ~gui.i_isuifig(host)
                % A mode left on owns the window's mouse callbacks.
                modes = {@datacursormode, @rotate3d, @zoom, @pan, @brush};
                for k = 1:numel(modes)
                    try
                        modes{k}(host, 'off');
                    catch
                        % mode not available for this figure; nothing to turn off
                    end
                end
            end

            if pending.Keep
                stack(end+1) = gui.myFigure.i_stashview(host, previous);
            else
                if ~isempty(previous) && isvalid(previous)
                    previous.i_release();
                end
                delete(host.Children);
            end
            setappdata(host, 'myFigureStack', stack);

            props = gui.myFigure.i_viewprops(host);
            for k = 1:numel(props)
                switch props{k}
                    case 'Name'
                        host.Name = '';
                    case 'Colormap'
                        host.Colormap = get(groot, 'DefaultFigureColormap');
                    case 'Pointer'
                        host.Pointer = 'arrow';
                    otherwise
                        host.(props{k}) = '';
                end
            end

            obj.FigHandle = host;
            obj.Adopted = true;
        end

        function i_suspend(obj)
            btns = [obj.tbv{~cellfun(@isempty, obj.tbv)}];
            btns = btns(isvalid(btns));
            obj.SuspendedButtons = btns(strcmp(get(btns, 'Visible'), 'on'));
            set(obj.SuspendedButtons, 'Visible', 'off');
        end

        function i_resume(obj)
            btns = obj.SuspendedButtons(isvalid(obj.SuspendedButtons));
            set(btns, 'Visible', 'on');
            obj.SuspendedButtons = gobjects(0);
        end

        function i_release(obj)
            % Let go of the window without closing it: the view's callbacks
            % die with its buttons, and with them this object, whose
            % destructor would otherwise close the window being reused.
            obj.FigHandle = [];
            btns = [obj.tbv{~cellfun(@isempty, obj.tbv)}];
            delete(btns(isvalid(btns)));
        end
    end

    methods (Static, Access = private)
        function pending = i_takehost(useuifig)
            pending = [];
            if ~isappdata(groot, 'myFigureHost'), return; end
            candidate = getappdata(groot, 'myFigureHost');
            % Claimed once, by the first myFigure DRAWFCN makes.
            gui.myFigure.i_clearhost();
            if pkg.i_isvalid(candidate.Fig) && isappdata(candidate.Fig, 'myFigureObj') ...
                    && gui.i_isuifig(candidate.Fig) == useuifig
                pending = candidate;
            end
        end

        function i_clearhost()
            if isappdata(groot, 'myFigureHost')
                rmappdata(groot, 'myFigureHost');
            end
        end

        function saved = i_stashview(host, previous)
            if ~isempty(previous) && isvalid(previous)
                previous.i_suspend();
            end
            % Hidden handle: HOST.CHILDREN, which the next view is built in
            % and which GOBACK clears, never includes it.
            panel = uipanel(host, 'Visible', 'off', 'BorderType', 'none', ...
                'Units', 'normalized', 'Position', [0, 0, 1, 1], ...
                'HandleVisibility', 'off', 'Tag', 'myFigureStash');
            kids = host.Children;
            hidden = gobjects(0);
            % Bottom first, so the panel keeps the stacking order.
            for k = numel(kids):-1:1
                kid = kids(k);
                switch kid.Type
                    case {'legend', 'colorbar'}
                        % They go with their axes.
                    case 'uicontextmenu'
                        % Never shown on its own; leave it with the figure.
                    case {'uitoolbar', 'uimenu'}
                        if strcmp(kid.Visible, 'on')
                            hidden(end+1) = kid; %#ok<AGROW>
                        end
                    otherwise
                        try
                            kid.Parent = panel;
                        catch
                            % Not allowed in a panel; hide it where it is.
                            if isprop(kid, 'Visible') && strcmp(kid.Visible, 'on')
                                hidden(end+1) = kid; %#ok<AGROW>
                            end
                        end
                end
            end
            % Out of HOST.CHILDREN too, or GOBACK would delete them along
            % with the view drawn over them.
            set(hidden, 'Visible', 'off', 'HandleVisibility', 'off');
            props = gui.myFigure.i_viewprops(host);
            figprops = cell2struct(get(host, props), props, 2);
            % Callbacks that act on "the plot in this window" read it.
            saved = struct('Obj', previous, 'Panel', panel, 'Hidden', hidden, ...
                'FigProps', figprops, 'CurrentAxes', host.CurrentAxes);
        end

        function props = i_viewprops(host)
            % ThemeChangedFcn is R2025a+.
            props = gui.myFigure.ViewProps;
            props = props(cellfun(@(x) isprop(host, x), props));
        end
    end
end
