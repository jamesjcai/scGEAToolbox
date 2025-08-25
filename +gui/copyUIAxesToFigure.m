function copyUIAxesToFigure(srcAx,ax2)
    % Create new figure and axes
    %fig2 = figure;
    %ax2 = axes(fig2);

    % Copy axes properties
    ax2.XLim        = srcAx.XLim;
    ax2.YLim        = srcAx.YLim;
    ax2.ZLim        = srcAx.ZLim;
    ax2.Title.String    = srcAx.Title.String;
    ax2.XLabel.String   = srcAx.XLabel.String;
    ax2.YLabel.String   = srcAx.YLabel.String;
    ax2.ZLabel.String   = srcAx.ZLabel.String;
    ax2.Box             = srcAx.Box;

    childObjs = srcAx.Children;

    % Copy contents
    for k = numel(childObjs):-1:1
        obj = childObjs(k);
        switch class(obj)
            case 'matlab.graphics.chart.primitive.Line'
                line(ax2, obj.XData, obj.YData, ...
                    'Color', obj.Color, 'LineStyle', obj.LineStyle, ...
                    'LineWidth', obj.LineWidth, 'Marker', obj.Marker, ...
                    'MarkerSize', obj.MarkerSize, 'DisplayName', obj.DisplayName);
            case 'matlab.graphics.chart.primitive.Scatter'
                scatter(ax2, obj.XData, obj.YData, ...
                    'CData', obj.CData, 'SizeData', obj.SizeData, ...
                    'Marker', obj.Marker, 'DisplayName', obj.DisplayName);
            case 'matlab.graphics.primitive.Image'
                image(ax2, 'XData', obj.XData, 'YData', obj.YData, 'CData', obj.CData);
            % Add other types if needed
        end
    end

    % Legends/colorbars
    if ~isempty(findobj(srcAx.Parent, 'Type', 'Legend'))
        legend(ax2);
    end
    copyColorbarComplete(srcAx, ax2);
    % copyColorbarBasic(srcAx, ax2);
    
    if ~isempty(findobj(srcAx.Parent, 'Type', 'ColorBar'))
        colorbar(ax2);    
    end
end


function copyColorbarComplete(srcAx, destAx)
% COPYCOLORBARCOMPLETE - Copy colorbar with all properties from source to destination axes
%
% Usage: copyColorbarComplete(srcAx, destAx)
%
% This function copies not just the colorbar, but all its properties including:
% - Label text and formatting
% - Tick marks and labels
% - Limits and direction
% - Location and size
% - Font properties
% - Colors and appearance

    % Find colorbar associated with source axes
    srcCbar = [];
    if isa(srcAx, 'matlab.ui.control.UIAxes')
        % For UIAxes, search in the parent figure
        allCbars = findobj(srcAx.Parent, 'Type', 'matlab.ui.control.ColorBar');
    else
        % For regular axes, search in the figure
        allCbars = findobj(srcAx.Parent, 'Type', 'ColorBar');
    end
    
    % Find the colorbar that belongs to our source axes
    for i = 1:length(allCbars)
        if isequal(allCbars(i).Axes, srcAx)
            srcCbar = allCbars(i);
            break;
        end
    end
    
    if isempty(srcCbar)
        % No colorbar found
        return;
    end
    
    % Create colorbar on destination axes
    destCbar = colorbar(destAx);
    
    % Copy all relevant properties
    try
        % Basic properties
        if isprop(srcCbar, 'Location') && isprop(destCbar, 'Location')
            destCbar.Location = srcCbar.Location;
        end
        
        if isprop(srcCbar, 'Position') && isprop(destCbar, 'Position')
            destCbar.Position = srcCbar.Position;
        end
        
        % Limits and direction
        if isprop(srcCbar, 'Limits') && isprop(destCbar, 'Limits')
            destCbar.Limits = srcCbar.Limits;
        end
        
        if isprop(srcCbar, 'Direction') && isprop(destCbar, 'Direction')
            destCbar.Direction = srcCbar.Direction;
        end
        
        % Ticks and tick labels
        if isprop(srcCbar, 'Ticks') && isprop(destCbar, 'Ticks')
            destCbar.Ticks = srcCbar.Ticks;
        end
        
        if isprop(srcCbar, 'TickLabels') && isprop(destCbar, 'TickLabels')
            destCbar.TickLabels = srcCbar.TickLabels;
        end
        
        if isprop(srcCbar, 'TickDirection') && isprop(destCbar, 'TickDirection')
            destCbar.TickDirection = srcCbar.TickDirection;
        end
        
        if isprop(srcCbar, 'TickLength') && isprop(destCbar, 'TickLength')
            destCbar.TickLength = srcCbar.TickLength;
        end
        
        % Label properties
        if isprop(srcCbar, 'Label') && isprop(destCbar, 'Label')
            % Copy label text
            destCbar.Label.String = srcCbar.Label.String;
            
            % Copy label font properties
            if isprop(srcCbar.Label, 'FontSize') && isprop(destCbar.Label, 'FontSize')
                destCbar.Label.FontSize = srcCbar.Label.FontSize;
            end
            
            if isprop(srcCbar.Label, 'FontWeight') && isprop(destCbar.Label, 'FontWeight')
                destCbar.Label.FontWeight = srcCbar.Label.FontWeight;
            end
            
            if isprop(srcCbar.Label, 'FontName') && isprop(destCbar.Label, 'FontName')
                destCbar.Label.FontName = srcCbar.Label.FontName;
            end
            
            if isprop(srcCbar.Label, 'FontAngle') && isprop(destCbar.Label, 'FontAngle')
                destCbar.Label.FontAngle = srcCbar.Label.FontAngle;
            end
            
            if isprop(srcCbar.Label, 'Color') && isprop(destCbar.Label, 'Color')
                destCbar.Label.Color = srcCbar.Label.Color;
            end
            
            if isprop(srcCbar.Label, 'Rotation') && isprop(destCbar.Label, 'Rotation')
                destCbar.Label.Rotation = srcCbar.Label.Rotation;
            end
        end
        
        % Font properties for the colorbar itself
        if isprop(srcCbar, 'FontSize') && isprop(destCbar, 'FontSize')
            destCbar.FontSize = srcCbar.FontSize;
        end
        
        if isprop(srcCbar, 'FontWeight') && isprop(destCbar, 'FontWeight')
            destCbar.FontWeight = srcCbar.FontWeight;
        end
        
        if isprop(srcCbar, 'FontName') && isprop(destCbar, 'FontName')
            destCbar.FontName = srcCbar.FontName;
        end
        
        if isprop(srcCbar, 'FontAngle') && isprop(destCbar, 'FontAngle')
            destCbar.FontAngle = srcCbar.FontAngle;
        end
        
        % Color properties
        if isprop(srcCbar, 'Color') && isprop(destCbar, 'Color')
            destCbar.Color = srcCbar.Color;
        end
        
        if isprop(srcCbar, 'Box') && isprop(destCbar, 'Box')
            destCbar.Box = srcCbar.Box;
        end
        
        % Additional appearance properties
        if isprop(srcCbar, 'LineWidth') && isprop(destCbar, 'LineWidth')
            destCbar.LineWidth = srcCbar.LineWidth;
        end
        
        if isprop(srcCbar, 'AxisLocation') && isprop(destCbar, 'AxisLocation')
            destCbar.AxisLocation = srcCbar.AxisLocation;
        end
        
        if isprop(srcCbar, 'AxisLocationMode') && isprop(destCbar, 'AxisLocationMode')
            destCbar.AxisLocationMode = srcCbar.AxisLocationMode;
        end
        
        % Ruler properties (for newer MATLAB versions)
        if isprop(srcCbar, 'Ruler') && isprop(destCbar, 'Ruler') && ~isempty(srcCbar.Ruler)
            try
                % Copy ruler properties if available
                if isprop(srcCbar.Ruler, 'TickLabelFormat') && isprop(destCbar.Ruler, 'TickLabelFormat')
                    destCbar.Ruler.TickLabelFormat = srcCbar.Ruler.TickLabelFormat;
                end
                
                if isprop(srcCbar.Ruler, 'Exponent') && isprop(destCbar.Ruler, 'Exponent')
                    destCbar.Ruler.Exponent = srcCbar.Ruler.Exponent;
                end
            catch
                % Ruler properties might not be available in all versions
            end
        end
        
    catch ME
        warning(ME.identifier, 'Some colorbar properties could not be copied: %s', ME.message);
    end
end

% Alternative simpler version for basic copying
function copyColorbarBasic(srcAx, destAx)
% COPYCOLORBARBASIC - Simple version that copies essential colorbar properties
    
    % Find source colorbar
    if isa(srcAx, 'matlab.ui.control.UIAxes')
        srcCbar = findobj(srcAx.Parent, 'Type', 'matlab.ui.control.ColorBar');
    else
        srcCbar = findobj(srcAx.Parent, 'Type', 'ColorBar');
    end
    
    if isempty(srcCbar)
        return;
    end
    
    % Handle multiple colorbars - find the one associated with srcAx
    for i = 1:length(srcCbar)
        if isequal(srcCbar(i).Axes, srcAx)
            srcCbar = srcCbar(i);
            break;
        end
    end
    
    if length(srcCbar) > 1
        srcCbar = srcCbar(1); % Take the first one if still multiple
    end
    
    % Create destination colorbar
    destCbar = colorbar(destAx);
    
    % Copy essential properties
    try
        destCbar.Label.String = srcCbar.Label.String;
        destCbar.Location = srcCbar.Location;
        destCbar.Ticks = srcCbar.Ticks;
        destCbar.TickLabels = srcCbar.TickLabels;
        destCbar.Limits = srcCbar.Limits;
        destCbar.FontSize = srcCbar.FontSize;
    catch
        warning('Some basic colorbar properties could not be copied');
    end
end