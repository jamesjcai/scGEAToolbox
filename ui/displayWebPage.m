function hContainer = displayWebPage(url, parent)
% displayWebPage - display a web-page URL in a Matlab figure or UI container
%
% Syntax: hContainer = displayWebPage(url, parent)
%
% displayWebPage() with no input parameters displays the UndocumentedMatlab.com
% homepage in a figure window.
%
% displayWebPage(url) displays the specified webpage in a figure window.
%
% displayWebPage(url, title) displays the specified webpage in a figure window
% that has the specified title, creating a new figure if no such figure is found.
% If the figure is found and already contains an embedded browser, its contents
% are updated with the new web-page. Otherwise, a new browser is added to the
% figure, and then loads the specified webpage.
%
% displayWebPage(url, hParent) displays the specified webpage in the specified
% container handle (figure, uipanel, uitab, etc.). If the specified handle is
% not valid, the webpage is displayed in the system browser.
%
% hContainer = displayWebPage(...) returns the browser's container handle
%
% Usage examples:
%    displayWebPage  % display UndocumentedMatlab.com homepage in a figure window
%    displayWebPage('google.com')  % display google.com in a figure window
%    displayWebPage('http://google.com', 'Browser')  % display in 'Browser' figure
%    displayWebPage('http://google.com', uipanel)    % display in specified panel
%    hContainer = displayWebPage('http://google.com');  % return the figure handle
%
% Notes:
%    In some cases when the specified webpage URL is invalid, the browser might
%    hang, causing excessive CPU load. Closing/deleting the browser container
%    (e.g. its figure window) will dispose the browser and restore CPU to normal.
%
% Additional information:
%    https://UndocumentedMatlab.com
%
% See also:
%    web

% Release history:
%    1.0  2022-10-01: initial version

% License to use and modify this code is granted freely to all interested,
% as long as the original author is referenced and attributed as such. 
% The original author maintains the right to be solely associated with this work.

% Programmed and Copyright by Yair M. Altman: altmany(at)gmail.com, UndocumentedMatlab.com

    try
        % Process missing optional input args
        if nargin < 1 || isempty(url)
            url = 'https://UndocumentedMatlab.com';
        end
        if nargin < 2 || isempty(parent)
            parent = 'Embedded browser figure window';
        end

        % Prepend 'http:' if no protocol was specified
        url = strtrim(char(url));
        if ~any(url==':'), url = ['http://' url]; end

        % If a container handle was specified, use it
        newFig = false;
        if ~ischar(parent) && ~isstring(parent)
            % Raise exception if the container handle is invalid
            if ~isvalid(parent), error('invalid container handle'); end

            % Get the browser's reference handle (if it exists) in the specified container handle
            try jBrowser = getappdata(parent,'jBrowser'); catch, jBrowser=[]; end

            % If container does not already contain embedded browser, create it
            if isempty(jBrowser), jBrowser = createBrowserIn(parent); end

            % Return the container's handle, if requested
            if nargout, hContainer = parent; end
        else
            % Check for existence of the figure with specified title
            hFig = findall(0, 'Name',parent, 'Tag','Browser figure', '-depth',1);
            try jBrowser = getappdata(hFig,'jBrowser'); catch, jBrowser=[]; end
            if isempty(hFig) || ~isvalid(hFig)
                % Create a new figure
                hFig = figure('Color','w', ... 'Units','norm', 'Pos',[0.3 0.2 0.4 0.5], ...
                    'Menubar','none', 'Toolbar','none', 'NumberTitle','off', ...
                    'Tag','Browser figure', 'Name',parent);
                newFig = true;

                % Add browser panel in the figure's main content pane
                jBrowser = createBrowserIn(hFig);

                % Add custom browser toolbar
                jAddressField = createBrowserToolbar(hFig, url);
            elseif isempty(jBrowser)  % figure exists but has no browser
                % Add browser panel in the figure's main content pane
                jBrowser = createBrowserIn(hFig);

                % Add custom browser toolbar
                jAddressField = createBrowserToolbar(hFig, url);
            else
                hFig = hFig(1);  % in case there are multiple matching figures
                %set(hFig,'Visible','on'); figure(hFig); %display & bring to focus

                % Get the browser's address field reference handle
                jAddressField = getappdata(hFig,'jAddressField');
            end

            % Return the figure's handle, if requested
            if nargout, hContainer = hFig; end

            % Update the URL address field (only in standalone figure, not panel)
            try jAddressField.setText(url); catch, end
        end

        % Load the specified URL in the embedded browser
        jBrowser.load(url);
    catch
        % Close any newly-created figure
        if newFig, delete(hFig), end

        % Open the URL in system browser
        web(url, '-browser');
        if nargout, hContainer = []; end
    end
end

% Create a new browser instance and place it in the specified container handle
function jBrowser = createBrowserIn(hParent)
    % Add maximized browser panel within hParent
    jBrowserPanel = javaObjectEDT(com.mathworks.mlwidgets.help.LightweightHelpPanel); %#ok<JAPIMATHWORKS>
    [jhBrowserPanel, hContainer] = javacomponent(jBrowserPanel, [], hParent); %#ok<JAVCM>
    set(hContainer, 'Units','norm', 'Position',[0,0,1,1]);

    % Store the browser reference handle for later use
    % Note: jBrowserPanel.setCurrentLocation(url) only displays URLs under
    % https://mathworks.com/help/, so we use jBrowser.load(url) instead
    jBrowser = jhBrowserPanel.getLightweightBrowser;
    setappdata(hParent,'jBrowser',jBrowser);

    % Set-up cleaner callback
    addlistener(hParent,'ObjectBeingDestroyed',@(h,e)cleanup(jBrowserPanel));
end

% Create a browser toolbar with an address-bar and an <open in browser> button
% Note: This toolbar is only created/displayed in standalone figure mode, not
% when the browser is embedded in an internal figure container (e.g. uipanel).
function jAddressField = createBrowserToolbar(hFig, url)
    % Add custom browser toolbar
    hToolbar = uitoolbar(hFig);
    drawnow  % required for jToolbar to be non-empty

    % Create the address box
    jAddressField = javaObjectEDT(javax.swing.JTextField(url));
    jAddressField.setEditable(false);
    %{
    jSize = java.awt.Dimension(400,25);
    jAddressField.setMaximumSize(jSize);
    jAddressField.setMinimumSize(jSize);
    jAddressField.setPreferredSize(jSize);
    jAddressField.setSize(jSize);
    %}
    setappdata(hFig,'jAddressField',jAddressField);

    % Add a narrow padding
    jPaddingPanel = javaObjectEDT(javax.swing.JPanel);
    jSize = java.awt.Dimension(3,25);
    jPaddingPanel.setMaximumSize(jSize);
    jPaddingPanel.setMinimumSize(jSize);
    jPaddingPanel.setPreferredSize(jSize);
    jPaddingPanel.setSize(jSize);

    % Create the simple push-button
    jOpenBrowser = javaObjectEDT(javax.swing.JButton('Open in browser'));
    jOpenBrowser.setToolTipText('Open this webpage in system browser');
    jhOpen = handle(jOpenBrowser, 'CallbackProperties');
    set(jhOpen, 'ActionPerformedCallback', @(h,e)web(char(jAddressField.getText)));

    % Append the filler and search-box to the toolbar
    jFiller = javax.swing.Box.createHorizontalGlue;  % javax.swing.Box$Filler
    jToolbar = hToolbar.JavaContainer.getComponentPeer; %#ok<JAVCT>
    jToolbar.add(jAddressField, jToolbar.getComponentCount);
    jToolbar.add(jFiller,       jToolbar.getComponentCount);
    jToolbar.add(jOpenBrowser,  jToolbar.getComponentCount);
    jToolbar.add(jPaddingPanel, jToolbar.getComponentCount);
    jToolbar.revalidate;
    jToolbar.repaint

    % Set-up cleaner callback to dispose browser resources upon container deletion
    addlistener(hFig,'ObjectBeingDestroyed',@(h,e)cleanup(jToolbar));
end

% Dispose the browser and its memory resources when Matlab container is deleted
function cleanup(jObject)
    try jObject.dispose; catch, end
end
