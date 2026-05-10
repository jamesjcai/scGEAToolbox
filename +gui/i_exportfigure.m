function i_exportfigure(h, filename)
%I_EXPORTFIGURE Export a figure or axes with fallbacks for constrained graphics.

if isa(h, 'matlab.ui.Figure')
    try
        exportapp(h, filename);
        return;
    catch
        % Fall through to alternate exporters.
    end
end

try
    exportgraphics(h, filename, 'BackgroundColor', 'white');
    return;
catch
    % Fall through to framebuffer capture.
end

try
    frame = getframe(h);
    imwrite(frame.cdata, filename);
    return;
catch ME
    error('scGEAToolbox:FigureExportFailed', ...
        'Unable to export figure to "%s": %s', filename, ME.message);
end
end
