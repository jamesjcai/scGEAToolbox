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

    properties
        FigHandle % Handle to the MATLAB figure
        AxHandle
    end

    properties (Access = private)
        tb
        tb2
        tbv = cell(12, 1) % Toolbar button handles
    end

    methods
        % Constructor
        function obj = myFigure(parentfig, forceclassic)
            if nargin<2, forceclassic = true; end
            if nargin<1, parentfig = []; end

            if ~forceclassic && gui.i_isuifig(parentfig)
                obj.FigHandle = uifigure('Name', '', ...
                    'Visible',"off","ToolBar","none");

                if ~isMATLABReleaseOlderThan('R2025a')
                    try
                        theme(obj.FigHandle, parentfig.Theme.BaseColorStyle);
                    catch
                    end
                end

                obj.tb = uitoolbar(obj.FigHandle);
                obj.AxHandle = uiaxes(obj.FigHandle);
                obj.AxHandle.Position = [50, 30, obj.FigHandle.Position(3)-100, ...
                    obj.FigHandle.Position(4)-60];

                obj.tbv{1} = gui.i_addbutton2fig(obj.tb, 'off', @gui.i_invertcolor, 'INVERT.gif', 'Invert Colors');
                obj.tbv{2} = gui.i_addbutton2fig(obj.tb, 'off', @gui.i_linksubplots, "keyframes-minus.jpg", "Link Subplots");
                obj.tbv{3} = gui.i_addbutton2fig(obj.tb, 'off', {@gui.i_setboxon, obj.FigHandle}, 'border-out.jpg', 'Box ON/OFF');
                obj.tbv{4} = gui.i_addbutton2fig(obj.tb, 'off', @gui.i_renametitle, "align-top-box.jpg", 'Add/Edit Title');
                obj.tbv{5} = gui.i_addbutton2fig(obj.tb, 'off', {@gui.i_pickcolor, false}, 'color-wheel.jpg', 'Pick a New Colormap...');
                obj.tbv{6} = gui.i_addbutton2fig(obj.tb, 'off', @gui.i_changefontsize, 'text-size.jpg', 'Change Font Size');
                obj.tbv{7} = gui.i_addbutton2fig(obj.tb, 'on', {@gui.i_savemainfig, 3}, "presentation.jpg", 'Save Figure to PowerPoint File...');
                obj.tbv{8} = gui.i_addbutton2fig(obj.tb, 'off', {@gui.i_savemainfig, 2, obj.FigHandle, obj.AxHandle}, "jpg-format.jpg", 'Save Figure as Graphic File...');
                obj.tbv{9} = gui.i_addbutton2fig(obj.tb, 'off', {@gui.i_savemainfig, 1, obj.FigHandle, obj.AxHandle}, "svg-format.jpg", 'Save Figure as SVG File...');
                obj.tbv{10} = gui.xui_3dcamera(obj.tb, '', false, obj.FigHandle, obj.AxHandle);
                obj.tbv{11} = gui.i_addbutton2fig(obj.tb, 'on', {@gui.i_resizewin, obj.FigHandle}, 'scale-frame-reduce.jpg', 'Resize Plot Window');
                obj.tbv{12} = gui.i_addbutton2fig(obj.tb, 'on', @obj.in_darkmode, 'demoIcon.gif', 'Light/Dark Mode');

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

                obj.AxHandle = axes('Parent', obj.FigHandle);
                obj.tb = findall(obj.FigHandle, 'Type', 'uitoolbar');
                if isempty(obj.tb)
                    obj.tb = uitoolbar(obj.FigHandle);
                end
                % Remove undesirable default tools if present
                delete(findall(obj.FigHandle, 'Tag', 'DataManager.Linking'));
                delete(findall(obj.FigHandle, 'Tag', 'Standard.OpenInspector'));
                obj.tbv{1} = gui.i_addbutton2fig(obj.tb, 'on', @gui.i_invertcolor, 'INVERT.gif', 'Invert Colors');
                obj.tbv{2} = gui.i_addbutton2fig(obj.tb, 'off', @gui.i_linksubplots, "keyframes-minus.jpg", "Link Subplots");
                obj.tbv{3} = gui.i_addbutton2fig(obj.tb, 'off', {@gui.i_setboxon, obj.FigHandle}, 'border-out.jpg', 'Box ON/OFF');
                obj.tbv{4} = gui.i_addbutton2fig(obj.tb, 'off', @gui.i_renametitle, "align-top-box.jpg", 'Add/Edit Title');
                obj.tbv{5} = gui.i_addbutton2fig(obj.tb, 'off', {@gui.i_pickcolor, false}, 'color-wheel.jpg', 'Pick a New Colormap...');
                obj.tbv{6} = gui.i_addbutton2fig(obj.tb, 'off', @gui.i_changefontsize, 'text-size.jpg', 'Change Font Size');
                obj.tbv{7} = gui.i_addbutton2fig(obj.tb, 'on', {@gui.i_savemainfig, 3}, "presentation.jpg", 'Save Figure to PowerPoint File...');
                obj.tbv{8} = gui.i_addbutton2fig(obj.tb, 'off', {@gui.i_savemainfig, 2, obj.FigHandle, obj.AxHandle}, "jpg-format.jpg", 'Save Figure as Graphic File...');
                obj.tbv{9} = gui.i_addbutton2fig(obj.tb, 'off', {@gui.i_savemainfig, 1, obj.FigHandle, obj.AxHandle}, "svg-format.jpg", 'Save Figure as SVG File...');
                obj.tbv{10} = gui.gui_3dcamera(obj.tb);
                obj.tbv{11} = gui.i_addbutton2fig(obj.tb, 'on', {@gui.i_resizewin, obj.FigHandle}, 'scale-frame-reduce.jpg', 'Resize Plot Window');
                obj.tbv{12} = gui.i_addbutton2fig(obj.tb, 'on', @obj.in_darkmode, 'demoIcon.gif', 'Light/Dark Mode');
            end

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
            centerto(obj, parentfig);
            obj.FigHandle.Visible = true;
        end

        function addCustomButton(obj, sepTag, callback, imgFil, tooltipTxt)
            if isempty(obj.tb2)
                obj.tb2 = uitoolbar(obj.FigHandle);
            end
            gui.i_addbutton2fig(obj.tb2, sepTag, callback, imgFil, tooltipTxt);
        end

        function setTitle(obj, titleStr)
            obj.FigHandle.Name = titleStr;
        end

        function closeFigure(obj)
            if isvalid(obj.FigHandle)
                close(obj.FigHandle);
            end
        end

        %% Destructor
        function delete(obj)
            closeFigure(obj);
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
    end
end
