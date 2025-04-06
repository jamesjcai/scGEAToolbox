classdef myFigure < handle
    properties
        FigHandle % Handle to the MATLAB figure        
        AxHandle
    end
    
    properties (Access = private)
        tb
        tb2
        tbv = cell(11,1)
    end

    methods
        % Constructor
        function obj = myFigure(parentfig)
            if nargin<1, parentfig = []; end

            if gui.i_isuifig(parentfig)
                % disp('making a uifigure');
                obj.FigHandle = uifigure('Name', '', ...
                    'Visible',"off");
                obj.tb = uitoolbar(obj.FigHandle);
                obj.AxHandle = uiaxes(obj.FigHandle);
                obj.AxHandle.Position = [50, 30, obj.FigHandle.Position(3)-100, ... 
                    obj.FigHandle.Position(4)-60]; % Fill most of the figure

                % pkg.i_addbutton2fig(obj.tb, 'on', [], [], '');
                obj.tbv{1} = pkg.i_addbutton2fig(obj.tb, 'off', @gui.i_invertcolor, 'INVERT.gif', 'Invert Colors');     
                obj.tbv{2} = pkg.i_addbutton2fig(obj.tb, 'off', @gui.i_linksubplots, "keyframes-minus.jpg", "Link Subplots");
                obj.tbv{3} = pkg.i_addbutton2fig(obj.tb, 'off', {@gui.i_setboxon, obj.FigHandle}, 'border-out.jpg', 'Box ON/OFF'); 
                obj.tbv{4} = pkg.i_addbutton2fig(obj.tb, 'off', @gui.i_renametitle, "align-top-box.jpg", 'Add/Edit Title');
                obj.tbv{5} = pkg.i_addbutton2fig(obj.tb, 'off', {@gui.i_pickcolor, false}, 'color-wheel.jpg', 'Pick a New Colormap...');
                obj.tbv{6} = pkg.i_addbutton2fig(obj.tb, 'off', @gui.i_changefontsize, 'text-size.jpg', 'Change Font Size');
                obj.tbv{7} = pkg.i_addbutton2fig(obj.tb, 'on', {@gui.i_savemainfig, 3}, "presentation.jpg", 'Save Figure to PowerPoint File...');
                obj.tbv{8} = pkg.i_addbutton2fig(obj.tb, 'off', {@gui.i_savemainfig, 2}, "jpg-format.jpg", 'Save Figure as Graphic File...');
                obj.tbv{9} = pkg.i_addbutton2fig(obj.tb, 'off', {@gui.i_savemainfig, 1}, "svg-format.jpg", 'Save Figure as SVG File...');
                obj.tbv{10} = gui.gui_3dcamera(obj.tb);
                obj.tbv{11} = pkg.i_addbutton2fig(obj.tb, 'on', {@gui.i_resizewin, obj.FigHandle}, 'scale-frame-reduce.jpg', 'Resize Plot Window');                
            else
                obj.FigHandle = figure('Name', '', ...
                    'NumberTitle', 'on', 'Visible',"off", ...
                    "DockControls", "off", 'ToolBar', 'figure', ...
                    'Position', [1, 1, 560, 420]);
                % obj.AxHandle = axes(obj.FigHandle);
                obj.tb = findall(obj.FigHandle, 'Type', 'uitoolbar');
                if isempty(obj.tb)
                    obj.tb = uitoolbar(obj.FigHandle);
                end
                linkTool = findall(obj.FigHandle, 'Tag', 'DataManager.Linking');
                if ~isempty(linkTool), delete(linkTool); end
                linkTool = findall(obj.FigHandle, 'Tag', 'Standard.OpenInspector');
                if ~isempty(linkTool), delete(linkTool); end
    
                % uipushtool(obj.tb, 'Separator', 'off');
                pkg.i_addbutton2fig(obj.tb, 'on', [], [], '');  
                obj.tbv{1} = pkg.i_addbutton2fig(obj.tb, 'on', @gui.i_invertcolor, 'INVERT.gif', 'Invert Colors');     
                obj.tbv{2} = pkg.i_addbutton2fig(obj.tb, 'off', @gui.i_linksubplots, "keyframes-minus.jpg", "Link Subplots");
                obj.tbv{3} = pkg.i_addbutton2fig(obj.tb, 'off', {@gui.i_setboxon, obj.FigHandle}, 'border-out.jpg', 'Box ON/OFF'); 
                obj.tbv{4} = pkg.i_addbutton2fig(obj.tb, 'off', @gui.i_renametitle, "align-top-box.jpg", 'Add/Edit Title');
                obj.tbv{5} = pkg.i_addbutton2fig(obj.tb, 'off', {@gui.i_pickcolor, false}, 'color-wheel.jpg', 'Pick a New Colormap...');
                obj.tbv{6} = pkg.i_addbutton2fig(obj.tb, 'off', @gui.i_changefontsize, 'text-size.jpg', 'Change Font Size');
                obj.tbv{7} = pkg.i_addbutton2fig(obj.tb, 'on', {@gui.i_savemainfig, 3}, "presentation.jpg", 'Save Figure to PowerPoint File...');
                obj.tbv{8} = pkg.i_addbutton2fig(obj.tb, 'off', {@gui.i_savemainfig, 2}, "jpg-format.jpg", 'Save Figure as Graphic File...');
                obj.tbv{9} = pkg.i_addbutton2fig(obj.tb, 'off', {@gui.i_savemainfig, 1}, "svg-format.jpg", 'Save Figure as SVG File...');
                obj.tbv{10} = gui.gui_3dcamera(obj.tb);
                obj.tbv{11} = pkg.i_addbutton2fig(obj.tb, 'on', {@gui.i_resizewin, obj.FigHandle}, 'scale-frame-reduce.jpg', 'Resize Plot Window');
            end
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
                    [px_new] = gui.i_getchildpos(parentfig, obj.FigHandle);
                    if ~isempty(px_new)
                        movegui(obj.FigHandle, px_new);
                    else
                        movegui(obj.FigHandle, 'center');
                    end
                else
                    movegui(obj.FigHandle, 'center');
                end
            catch
                movegui(obj.FigHandle, 'center');
            end        
        end

        function show(obj, parentfig)
            if nargin<2, parentfig = []; end
            centerto(obj, parentfig);
            obj.FigHandle.Visible = true;
        end
        
        function addCustomButton(obj, sepTag, callback, imgFil, buttonLabel)
            if isempty(obj.tb2)
                obj.tb2 = uitoolbar(obj.FigHandle);
            end
            pkg.i_addbutton2fig(obj.tb2, sepTag, callback, imgFil, buttonLabel); 
        end

        function setTitle(obj, titleStr)
            obj.FigHandle.Name = titleStr;
        end
        
        function closeFigure(obj)
            close(obj.FigHandle);
        end
    end
end
