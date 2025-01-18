classdef myFigure < handle
    properties
        FigureHandle % Handle to the MATLAB figure        
    end
    
    properties (Access = private)
        tb
        tb2
        tbv = cell(11,1)
    end

    methods
        % Constructor
        function obj = myFigure()
            obj.FigureHandle = figure('Name', '', ...
                'NumberTitle', 'on', 'Visible',"off", ...
                "DockControls", "off", 'ToolBar', 'figure', ...
                'Position', [1, 1, 560, 420]);
            obj.tb = findall(obj.FigureHandle, 'Type', 'uitoolbar');
            if isempty(obj.tb)
                obj.tb = uitoolbar(obj.FigureHandle);
            end

            linkTool = findall(obj.FigureHandle, 'Tag', 'DataManager.Linking');
            if ~isempty(linkTool), delete(linkTool); end
            linkTool = findall(obj.FigureHandle, 'Tag', 'Standard.OpenInspector');
            if ~isempty(linkTool), delete(linkTool); end

            % uipushtool(obj.tb, 'Separator', 'off');
            pkg.i_addbutton2fig(obj.tb, 'on', [], [], '');  
            obj.tbv{1} = pkg.i_addbutton2fig(obj.tb, 'on', @gui.i_invertcolor, 'INVERT.gif', 'Invert colors');     
            obj.tbv{2} = pkg.i_addbutton2fig(obj.tb, 'off', @gui.i_linksubplots, "keyframes-minus.jpg", "Link subplots");
            obj.tbv{3} = pkg.i_addbutton2fig(obj.tb, 'off', @gui.i_setboxon, 'border-out.jpg', 'Box on/off'); 
            obj.tbv{4} = pkg.i_addbutton2fig(obj.tb, 'off', @gui.i_renametitle, "align-top-box.jpg", 'Change Plot Title');
            obj.tbv{5} = pkg.i_addbutton2fig(obj.tb, 'off', {@gui.i_pickmonocolor, true}, 'color-wheel.jpg', 'Pick new color map...');
            obj.tbv{6} = pkg.i_addbutton2fig(obj.tb, 'off', @gui.i_changefontsize, 'text-size.jpg', 'ChangeFontSize');
            obj.tbv{7} = pkg.i_addbutton2fig(obj.tb, 'on', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
            obj.tbv{8} = pkg.i_addbutton2fig(obj.tb, 'off', {@gui.i_savemainfig, 2}, "File-Jpg--Streamline-Core-Gradient.png", 'Save Figure as Graphic File...');
            obj.tbv{9} = pkg.i_addbutton2fig(obj.tb, 'off', {@gui.i_savemainfig, 1}, "svg.png", 'Save Figure as SVG File...');
            obj.tbv{10} = gui.gui_3dcamera(obj.tb, 'AllCells');
            obj.tbv{11} = pkg.i_addbutton2fig(obj.tb, 'on', {@gui.i_resizewin, obj.FigureHandle}, 'scale-frame-enlarge.jpg', 'Resize Plot Window');
        end

        % function ptvkeep(obj, flag)
        %     delete(obj.tbv{~flag});
        % end

        function ptvshow(obj, flag)
            set([obj.tbv{flag}],'Visible','on');
            set([obj.tbv{~flag}],'Visible','off');
        end

        function ptvenable(obj, flag)
            set([obj.tbv{flag}],'Enable','on');
            set([obj.tbv{~flag}],'Enable','off');
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

        function show(obj, parentfig)
            if nargin<2, parentfig = []; end
            centerto(obj, parentfig);
            obj.FigureHandle.Visible = true;
        end
        
        function addCustomButton(obj, sepTag, callback, imgFil, buttonLabel)
            if isempty(obj.tb2)
                obj.tb2 = uitoolbar(obj.FigureHandle);
            end
            pkg.i_addbutton2fig(obj.tb2, sepTag, callback, imgFil, buttonLabel); 
        end

        function setTitle(obj, titleStr)
            obj.FigureHandle.Name = titleStr;
        end
        
        function closeFigure(obj)
            close(obj.FigureHandle);
        end
    end
end
