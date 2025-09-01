function i_savemainfig(src, ~, tag, fig, ax)

    % Initialize the figure and axes handles
    if nargin < 4, fig = []; end
    if nargin < 5, ax = []; end

    if ~isempty(fig)
        parentfig = fig;
    else
        parentfig = gui.gui_getfigsce(src);    
    end
    if ~isempty(ax)
        axesHandles = ax;
    else
        axesHandles = findall(parentfig, 'Type', 'axes');
    end

    if isempty(axesHandles)
        gui.myHelpdlg(parentfig, ...
            'No figures available in the current window. Unable to save figure.');
        return;
    end

    switch tag
        case 1 % PDF only
            defFilter = {'*.pdf'};
            [filename, filepath, filterIdx] = uiputfile(defFilter, 'Save as PDF');
            if isequal(filename,0), return; end
            fullfileOut = enforceExtension(fullfile(filepath, filename), defFilter{filterIdx});

            if gui.i_isuifig(parentfig) || ~isMATLABReleaseOlderThan('R2025a')
                exportWhite(parentfig, fullfileOut);
            else
                saveas(parentfig, fullfileOut, 'pdf');
            end

        case 2 % Multiple formats
            filters = {'*.jpg'; '*.png'; '*.tif'; '*.pdf'; '*.eps'};
            [filename, filepath, filterIdx] = uiputfile(filters, 'Save Figure');
            if isequal(filename,0), return; end
            fullfileOut = enforceExtension(fullfile(filepath, filename), filters{filterIdx});

            if gui.i_isuifig(parentfig) || ~isMATLABReleaseOlderThan('R2025a')
                exportWhite(parentfig, fullfileOut);
            else
                exportgraphics(parentfig, fullfileOut);
            end

        case 3 % Export to PPTX
            gui.i_export2pptx({parentfig}, {''}, parentfig);
    end
end


function outFile = enforceExtension(filepath, filterPattern)
    % Ensure the filename has the extension from the filter
    [~,~,ext] = fileparts(filepath);
    if isempty(ext)
        % filterPattern looks like '*.pdf' → get '.pdf'
        [~,~,ext] = fileparts(filterPattern);
        outFile = [filepath ext];
    else
        outFile = filepath;
    end
end

function exportWhite(h, filename)
    % exportWhite Export figure or axes with a temporary white background
    %
    %   exportWhite(h, filename)
    %
    %   h can be a matlab.ui.Figure (uifigure) or Axes handle.

    state = struct('h',{},'prop',{},'val',{});
    push = @(obj,prop,val) addState(obj,prop,val);        

    function addState(obj, propName, ~)
        try
            oldVal = get(obj, propName);
        catch
            return
        end
        state(end+1).h = obj; %#ok<AGROW>
        state(end).prop = propName;
        state(end).val = oldVal;

        try
            if isnumeric(oldVal) && size(oldVal,2) == 3 && size(oldVal,1) > 1
                set(obj, propName, repmat([1 1 1], size(oldVal,1), 1));
            else
                set(obj, propName, [1 1 1]);
            end
        catch
            % ignore components that reject the assignment
        end
    end    

    cleaner = onCleanup(@() restoreState(state));

    if isa(h, 'matlab.ui.Figure')
        theme(h,"light");
        allObjs = [h; findall(h)];
    for k = 1:numel(allObjs)
        obj = allObjs(k);
        if isprop(obj,'BackgroundColor')
            % Avoid colorbar (text color confusion) – it doesn't use BackgroundColor anyway.
            % Most UI components accept BackgroundColor=[1 1 1]; restore later.
            push(obj,'BackgroundColor',[]);
        end
    end        

        % --- Case 1: UIFigure (App Designer)
        oldColor = h.Color;
        cleanupObj = onCleanup(@() set(h, 'Color', oldColor)); 

        h.Color = [1 1 1];
        exportapp(h, filename);
    elseif isa(h, 'matlab.graphics.axis.Axes') || isa(h, 'matlab.ui.Figure')
        % --- Case 2: Regular Figure or Axes
        oldColor = h.Color;
        cleanupObj = onCleanup(@() set(h, 'Color', oldColor)); 

        h.Color = [1 1 1];
        exportgraphics(h, filename, 'BackgroundColor','white');

    else
        error('Unsupported handle type: must be a Figure or Axes.');
    end
end



%{
function exportWhite_new(h, filename)
% exportWhite Export a figure/axes/app with a temporarily white background.
%   exportWhite(h, filename)
%   - h can be a uifigure/figure/axes (UIAxes or regular Axes).
%   - All relevant containers & axes are set to white, exported, then restored.

    % Find everything under h (including h)
    allObjs = [h; findall(h)];
    state = struct('h',{},'prop',{},'val',{});
    push = @(obj,prop,val) addState(obj,prop,val);
    % Helper to push old value & set to white (handles Nx3 table BG too)
    function addState(obj, propName, ~)
        try
            oldVal = get(obj, propName);
        catch
            return
        end
        state(end+1).h = obj; %#ok<AGROW>
        state(end).prop = propName;
        state(end).val = oldVal;

        try
            if isnumeric(oldVal) && size(oldVal,2) == 3 && size(oldVal,1) > 1
                set(obj, propName, repmat([1 1 1], size(oldVal,1), 1));
            else
                set(obj, propName, [1 1 1]);
            end
        catch
            % ignore components that reject the assignment
        end
    end    

    % Helper to push old value & set to white (handles Nx3 table BG too)
    function addState(obj, propName, ~)
        try
            oldVal = get(obj, propName);
        catch
            return
        end
        state(end+1).h = obj; %#ok<AGROW>
        state(end).prop = propName;
        state(end).val = oldVal;

        try
            if isnumeric(oldVal) && size(oldVal,2) == 3 && size(oldVal,1) > 1
                set(obj, propName, repmat([1 1 1], size(oldVal,1), 1));
            else
                set(obj, propName, [1 1 1]);
            end
        catch
            % ignore components that reject the assignment
        end
    end

    % Identify types safely
    function tf = isAxesLike(obj)
        t = get(obj,'Type');
        tf = any(strcmp(t, {'axes','uiaxes','polaraxes','geoaxes'}));
    end

    % 1) Figure/app background
    for k = 1:numel(allObjs)
        obj = allObjs(k);
        if isgraphics(obj,'figure')
            if isprop(obj,'Color'), push(obj,'Color',[]); end
        end
    end

    % 2) Axes-like (Axes, UIAxes, polaraxes, geoaxes) → Color to white
    for k = 1:numel(allObjs)
        obj = allObjs(k);
        if isAxesLike(obj) && isprop(obj,'Color')
            push(obj,'Color',[]);
        end
    end

    % 3) Legends (background)
    for k = 1:numel(allObjs)
        obj = allObjs(k);
        if isgraphics(obj,'legend') && isprop(obj,'Color')
            push(obj,'Color',[]);
        end
    end

    % 4) Layout/containers & UI components with BackgroundColor
    %    (uitab, uitabgroup, uipanel, uigridlayout, and most uicontrols)
    for k = 1:numel(allObjs)
        obj = allObjs(k);
        if isprop(obj,'BackgroundColor')
            % Avoid colorbar (text color confusion) – it doesn't use BackgroundColor anyway.
            % Most UI components accept BackgroundColor=[1 1 1]; restore later.
            push(obj,'BackgroundColor',[]);
        end
    end

    % Ensure restoration no matter what
    cleaner = onCleanup(@() restoreState(state));

    % Export using the appropriate function
    if isa(h,'matlab.ui.Figure')
        % uifigure → exportapp (fallback to exportgraphics if unavailable)
        try
            exportapp(h, filename);
        catch
            exportgraphics(h, filename, 'BackgroundColor','white');
        end
    elseif isgraphics(h,'axes')
        exportgraphics(h, filename, 'BackgroundColor','white');
    else
        % classic figure (or other exportgraphics-supported container)
        exportgraphics(h, filename, 'BackgroundColor','white');
    end
end
%}

function restoreState(state)
% Restore all saved properties in reverse order
    for k = numel(state):-1:1
        s = state(k);
        if isvalid(s.h)
            try
                set(s.h, s.prop, s.val);
            catch
                % ignore if component disappeared or prop changed
            end
        end
    end
end
