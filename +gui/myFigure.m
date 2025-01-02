classdef myFigure < handle
    properties
        FigureHandle % Handle to the MATLAB figure
    end
    
    methods
        % Constructor
        function obj = myFigure()
            obj.FigureHandle = figure('Name', 'My Custom Figure', 'NumberTitle', 'off', 'Visible',"off");
            if isempty(findall(obj.FigureHandle, 'Type', 'uitoolbar'))
                tb = uitoolbar(obj.FigureHandle);
            else
                tb = findall(obj.FigureHandle, 'Type', 'uitoolbar');
            end
            pkg.i_addbutton2fig(tb, 'off', @gui.i_invertcolor, 'plotpicker-comet.gif', 'Invert colors');     
            pkg.i_addbutton2fig(tb, 'off', @gui.i_linksubplots, "plottypectl-rlocusplot.gif", "Link subplots");
            pkg.i_addbutton2fig(tb, 'off', @gui.i_setboxon, 'RectGate.gif', 'Box on/off'); 
            pkg.i_addbutton2fig(tb, 'off', @gui.i_renametitle, "icon-mat-touch-app-10.gif", 'Change Plot Title');
            pkg.i_addbutton2fig(tb, 'off', {@gui.i_pickmonocolor, true}, 'plotpicker-compass.gif', 'Pick new color map...');
            pkg.i_addbutton2fig(tb, 'off', @gui.i_changefontsize, 'noun_font_size_591141.gif', 'ChangeFontSize');
            pkg.i_addbutton2fig(tb, 'on', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
            pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 2}, "svg.gif", 'Save Figure as Graphic File...');
            pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 1}, "svg.gif", 'Save Figure as SVG File...');
            pkg.i_addbutton2fig(tb, 'on', {@gui.i_resizewin, obj.FigureHandle}, 'HDF_pointx.gif', 'Resize Plot Window');
            gui.gui_3dcamera(tb, 'AllCells');
            
        end

        function centerto(obj, parentfig)
            try
                if ~isempty(parentfig) && isa(parentfig,'matlab.ui.Figure') 
                    [px_new] = gui.i_getchildpos(parentfig, obj.FigureHandle);
                    if ~isempty(px_new)
                        movegui(obj.FigureHandle, px_new);
                    else
                        movegui(obj.FigureHandle, 'center');
                    end
                else
                    movegui(obj.FigureHandle, 'center');
                end
            catch
                movegui(obj.FigureHandle, 'center');
            end        
        end

        function show(obj)
            obj.FigureHandle.Visible = true;
        end
        
        % Method to add a custom toolbar button
        function addCustomButton(obj, buttonLabel, callback)
            % Ensure the figure has a toolbar
            if isempty(findall(obj.FigureHandle, 'Type', 'uitoolbar'))
                tb = uitoolbar(obj.FigureHandle);
            else
                tb = findall(obj.FigureHandle, 'Type', 'uitoolbar');
            end
            
            % Create a push tool
            pt = uipushtool(tb, 'TooltipString', buttonLabel);
            
            % Set the icon (optional)
            % pt.CData = imread('icon.png'); % Ensure the icon is a 16x16 RGB image
            pt.CData = rand(16, 16, 3);
            % Set the callback function
            pt.ClickedCallback = callback;
        end
        
        % Method to set the figure title
        function setTitle(obj, titleStr)
            obj.FigureHandle.Name = titleStr;
        end
        
        % Method to close the figure
        function closeFigure(obj)
            close(obj.FigureHandle);
        end
    end
end
