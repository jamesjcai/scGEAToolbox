function i_singlegraph(G1, figname, parentfig)

if nargin < 3, parentfig = []; end
if nargin < 2, figname = ''; end
if nargin < 1
    G1 = WattsStrogatz(100, 5, 0.15);
    G1.Nodes.Name = string((1:100)');
    G1.Edges.Weight = rand(size(G1.Edges, 1), 1) * 2;
end
import gui.*
import ten.*

%%

mfolder = fileparts(mfilename('fullpath'));

load(fullfile(mfolder, ...
    '..', 'resources', 'TFome', 'tfome_tfgenes.mat'), 'tfgenes');


w = 3;
l = 1;


hx = gui.myFigure;
hFig = hx.FigureHandle;
h1 = axes(hFig);
xy = [];
[p1] = drawnetwork(G1, h1);


%hx.ptvenable(logical([1 1 1 1 1 0 1 1 1 1 0]));

hx.addCustomButton('off', @ChangeFontSize, 'noun_font_size_591141.gif', 'Change Font Size of Nodes');
hx.addCustomButton('off', @ChangeWeight, 'weight_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Change Width of Edges');
hx.addCustomButton('off', @ChangeLayout, 'group_work_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Change Network Layout');
hx.addCustomButton('off', @ChangeDirected, 'turn_sharp_right_17dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Change to Undirected Network');
hx.addCustomButton('off', @AnimateCutoff, 'movie.jpg', 'Animate Change and Select a Cutoff for Linked Edges');
hx.addCustomButton('off', @ChangeCutoff, 'carpenter_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Change Cutoff to Trim Network');
hx.addCustomButton('off', @SaveAdj, 'floppy-disk-arrow-in.jpg', 'Export & Save Data');
hx.addCustomButton('on',  @in_RefreshAll, "refresh.jpg", "Refresh View");
hx.addCustomButton( 'on', @in_networkvis_linear, "linear.jpg", "Stright Network Vis");
hx.addCustomButton( 'off', @in_networkvis_curvy, "curve-array.jpg", "Curvy Network Vis");

title(h1,figname);
hx.show(parentfig);
% gui.gui_showrefinfo('Network Legend');

oldG1=[];
axistrig = true;
% dataengated = false;
oldidx = 0;



    function [xy] = getxy
        xy = [p1.XData' p1.YData'];
    end

    function in_networkvis_curvy(~, ~)
        fw=gui.gui_waitbar;
        gui.i_networkvis(G1, getxy, true, p1.NodeFontSize, hFig);
        gui.gui_waitbar(fw);
    end


    function in_networkvis_linear(~, ~)
        fw=gui.gui_waitbar;
        %figure('SizeChangedFcn', @sbar);
        h = gui.myFigure;        
        %assignin('base',"XYCoords",[p1.XData' p1.YData']);
        
        [x, y] = gplot(G1.adjacency, [p1.XData' p1.YData']);
        %gplot(G1.adjacency, [p1.XData' p1.YData'], '-k');
        plot(x, y,'k-');
        hold on
        % plot(p(:,1), p(:,2),'k^');
        % quiver(p(4,1),p(4,2),p(4,1)+0.1,p(4,2)+0.1,'r')
        % quiver(p(1,1),p(1,2),p(1,1)+0.1,p(1,2)+0.1,'r')         
        % customeMarker_1(p(:,1), p(:,2), theta_deg, h.FigureHandle);
        customeMarker_2(x, y, h.FigureHandle);

        % scatter(p1.XData', p1.YData', 300, ...
        %     'MarkerEdgeColor','k', ...
        %     'MarkerFaceColor',[.8 .8 .8]);

        textOpts.FontSize = p1.NodeFontSize;
        textOpts.HorizontalAlignment = 'center';
        textOpts.VerticalAlignment = 'middle';
        textOpts.FontWeight = 'normal';
        
        tz = cell(length(G1.Nodes.Name),1);
        for k = 1:length(G1.Nodes.Name)            
            [wx] = measureText(G1.Nodes.Name{k}, textOpts);
            tz{k} = text(p1.XData(k)-floor(wx/2), p1.YData(k), ...
                G1.Nodes.Name{k},'FontSize',textOpts.FontSize,...
                'BackgroundColor','w', ...
                'FontWeight','normal', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle');             
        end
        set(gca, 'XTick', [], 'YTick', []);
        axis off
        gui.gui_waitbar(fw);
        h.show(hFig);
        set(gcf, 'Color', 'white');
    end

    % function in_NetworkVis(~, ~)
    %     net = gui.networkvisualizer(G1.adjacency);
    %     net.setNodeLabels(G1.Nodes.Name);
    %     net.setNodeSizes('auto');
    %     net.X = p1.XData';
    %     net.Y = p1.YData';
    %     figure;
    %     plot(net);
    % end

    function in_RefreshAll(~, ~)
        if ~isempty(oldG1)
            G1 = oldG1;
        end
        [p1] = drawnetwork(G1, h1);
        title(h1,figname);
    end

    function SaveAdj(~, ~)
        if ~(ismcc || isdeployed)
            answer = questdlg('Export & save network to:', '', ...
                'Workspace', 'File', 'Workspace');
        else
            a = questdlg('Export & save network to file?');
            if strcmp(a, 'Yes')
                answer = 'File';
            else
                return;
            end
        end

        switch answer
            case 'Workspace'
                labels = {'Save adjacency matrix A1 to variable named:', ...
                    'Save graph G1 to variable named:', ...
                    'Save genelist g1 to variable named:'};
                A1 = full(adjacency(G1, 'weighted'));
                g1 = string(G1.Nodes.Name);
                vars = {'A', 'G', 'g'}; ...
                    values = {A1, G1, g1};
                msgfig = export2wsdlg(labels, vars, values);
                uiwait(msgfig);
            case 'File'
                [file, path] = uiputfile({'*.txt'; '*.*'}, 'Save as');
                if isequal(file, 0) || isequal(path, 0)
                    return;
                else
                    filename = fullfile(path, file);
                    pkg.m2cyto(G1, filename);
                end
            otherwise
                return;
        end
    end

    function ChangeFontSize(~, ~)
        i_changefontsize(p1);
        %i_changefontsize(p2);
        function i_changefontsize(p)
            if p.NodeFontSize >= 20
                p.NodeFontSize = 7;
            else
                p.NodeFontSize = p.NodeFontSize + 1;
            end
        end            
    end

    function ChangeBox(~, ~)
        if h1.Box
            box(h1, 'off');
            if axistrig
                axis(h1, 'off');
            else
                axis(h1, 'on');
            end
            axistrig = ~axistrig;
            
        else
            box(h1, 'on');
            axis(h1, 'on');
        end
    end

    function ChangeWeight(~, ~)
        %a=3:10;
        %w=a(randi(length(a),1));
        w = w + 1;
        if w > 10, w = 2; end
        i_changeweight(p1, w);
        %i_changeweight(p2,G2,w);
        function i_changeweight(p, b)
            %G.Edges.LWidths = abs(b*G.Edges.Weight/max(abs(G.Edges.Weight)));
            %p.LineWidth = G.Edges.LWidths;
            p.LineWidth = abs(b*p.LineWidth/max(abs(p.LineWidth)));
        end
    end

    function ChangeLayout(~, ~)
        a = ["auto", "layered", "subspace", "force", "circle", "reset"];
        l = l + 1;
        if l > 5, l = 1; end
        switch a(l)
            case "force"
                p1.layout(a(l), 'Iterations', 250);
            case "reset"
                [p1] = drawnetwork(G1, h1);
            otherwise
                p1.layout(a(l));
        end
    end

    function ChangeDirected(~, ~)
        if isempty(oldG1)
            oldG1 = G1;
        end
        [p1, G1] = i_changedirected(p1, G1, h1);

        function [p, G] = i_changedirected(p, G, h)
            x = p.XData;
            y = p.YData;
            if isa(G, 'digraph')
                A = adjacency(G, 'weighted');
                G = graph(0.5*(A + A.'), G.Nodes.Name);
                % p=plot(h,G);
                [p] = drawnetwork(G, h);
            end
            p.XData = x;
            p.YData = y;
        end
    end
                            
    function ChangeCutoff(~, ~)
        list = {'0.00 (show all edges)', ...
            '0.30', '0.35', '0.40', '0.45', ...
            '0.50', '0.55', '0.60', ...
            '0.65', '0.70', '0.75', '0.80', '0.85', ...
            '0.90', '0.95 (show 5% of edges)'};
        [indx, tf] = listdlg('ListString', list, ...
            'SelectionMode', 'single', ...
            'ListSize', [220, 300]);
        if tf
            if indx == 1
                cutoff = 0;
            elseif indx == length(list)
                cutoff = 0.95;
            else
                cutoff = str2double(list(indx));
            end
            answer = questdlg("Keep original network?","");
            switch answer
                case 'Yes'
            [p1] = i_replotg(p1, G1, h1, cutoff);
                case 'No'
            [p1, G1] = i_replotg(p1, G1, h1, cutoff);
                otherwise
                    return;
            end
            %[p2]=i_replotg(p2,G2,h2,cutoff);
        end
    end

    function [p] = drawnetwork(G, h)
        p = plot(h, G, 'ButtonDownFcn', @startDragFcn);
        box off
        layout(p,'force');
        %         if isa(G,'digraph')
        %             G.Nodes.NodeColors = outdegree(G)-indegree(G);
        %         else
        %             G.Nodes.NodeColors = degree(G);
        %         end
        %         p.NodeCData = G.Nodes.NodeColors;

        cc = repmat([0, 0.4470, 0.7410], G.numedges, 1);
        cc(G.Edges.Weight < 0, :) = repmat([0.8500, 0.3250, 0.0980], ...
            sum(G.Edges.Weight < 0), 1);
        p.EdgeColor = cc;
        %       p.EdgeCData=ones(G.numedges,1);
        %       p.EdgeCData(G.Edges.Weight<0)=2;

        ix = ismember(string(upper(G.Nodes.Name)), tfgenes);
        if any(ix)
            cc = repmat([0, 0, 0], G.numnodes, 1);
            cc(ix, :) = repmat([1, 0, 0], sum(ix), 1);
            p.NodeLabelColor = cc;
        end

        %p.NodeFontSize = 2 * p.NodeFontSize;

        %title(h,sprintf('%d nodes',G.numnodes));
        % https://www.mathworks.com/matlabcentral/answers/296070-change-label-font-in-graph-plots
        %{
        nl = p.NodeLabel;
        p.NodeLabel = '';
        xd = get(p, 'XData');
        yd = get(p, 'YData');
        text(xd, yd, nl, 'FontSize',p.NodeFontSize,...
            'FontWeight','bold',...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','middle',...
            'BackgroundColor','w','Margin',0.1);
        %}
        if ~isempty(G.Edges.Weight)
            G.Edges.LWidths = abs(w*G.Edges.Weight/max(G.Edges.Weight));
            p.LineWidth = G.Edges.LWidths;
        end
        xy = [p.XData' p.YData'];
    end

    function AnimateCutoff(~, ~)
        listc = 0.05:0.05:0.95;
        f = waitbar(0, 'Cutoff = 0.05', 'Name', 'Edge Pruning...', ...
            'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
        setappdata(f, 'canceling', 0);
    
        m = length(listc);
        for k = 1:m
            if getappdata(f, 'canceling')
                break
            end
    
            cutoff = listc(k);
            waitbar(k/m, f, sprintf('Cutoff = %g', cutoff));
            try
                p1 = i_replotg(p1, G1, h1, cutoff);
            catch ME
                disp(ME.message);
            end
            %p2=i_replotg(p2,G2,h2,cutoff);
            pause(1);
        end
        %close(f)
        delete(f)
    end

    function [p, G] = i_replotg(p, G, h, cutoff)
        a = h.Title.String;
        x = p.XData;
        y = p.YData;
        A = adjacency(G, 'weighted');
        A = ten.e_filtadjc(A, cutoff);
        if issymmetric(A)
            G = graph(A, G.Nodes.Name);
        else
            G = digraph(A, G.Nodes.Name);
        end
        % p=plot(h,G);
        [p] = drawnetwork(G, h);
        p.XData = x;
        p.YData = y;
        % h=gca;
        title(a)
    end

    %{
    % Callback to initiate dragging
    function startDragFcn(hObj, ~)
        % Get data for the clicked point
        fig = ancestor(hObj, 'figure');
        set(fig, 'WindowButtonMotionFcn', {@draggingFcn, hObj});
        set(fig, 'WindowButtonUpFcn', @stopDragFcn);
    end
    
    % Function to drag the point
    function draggingFcn(~, ~, hObj)
        % Current cursor position in data coordinates
        cp = get(gca, 'CurrentPoint');
        % Update the y-data of the nearest point
        yData = get(hObj, 'YData');
        xData = get(hObj, 'XData');
        [~, idx] = min(abs(xData - cp(1,1))); % Find closest x to mouse
        yData(idx) = cp(1,2); % Update y value
        set(hObj, 'YData', yData);
    end
    
    % Function to stop dragging
    function stopDragFcn(~, ~)
        fig = gcbf;
        set(fig, 'WindowButtonMotionFcn', '');
        set(fig, 'WindowButtonUpFcn', '');
    end
    %}

    % Callback to initiate dragging
    function startDragFcn(hObj, ~)
        % Get data for the clicked point
        % fig = ancestor(hObj, 'figure');
        set(hFig, 'WindowButtonMotionFcn', {@draggingFcn, hObj});
        set(hFig, 'WindowButtonUpFcn', @stopDragFcn);
    end
    
    % Function to drag the point
    function draggingFcn(~, ~, hObj)

            % Current cursor position in data coordinates
            cp = get(gca, 'CurrentPoint');
            % Update the y-data of the nearest point
            yData = get(hObj, 'YData');
            xData = get(hObj, 'XData');
    %       [~, idx] = min(abs(xData - cp(1,1))); % Find closest x to mouse
    %       yData(idx) = cp(1,2); % Update y value
    %       set(hObj, 'YData', yData);
            if oldidx == 0 % ~dataengated
                idx = dsearchn([xData' yData'], [cp(1,1) cp(1,2)]);
                oldidx = idx;
                % dataengated = true;
            else
                idx = oldidx;
            end
            xData(idx) = cp(1,1);
            yData(idx) = cp(1,2);
            set(hObj, 'XData', xData); % Update y value
            set(hObj, 'YData', yData);
    end
    
    % Function to stop dragging
    function stopDragFcn(~, ~)
        % fig = gcbf;
        set(hFig, 'WindowButtonMotionFcn', '');
        set(hFig, 'WindowButtonUpFcn', '');
        % dataengated = false;
        oldidx = 0;
    end


end



function h = WattsStrogatz(N, K, beta)
        % H = WattsStrogatz(N,K,beta) returns a Watts-Strogatz model graph with N
        % nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
        %
        % beta = 0 is a ring lattice, and beta = 1 is a random graph.

        % Connect each node to its K next and previous neighbors. This constructs
        % indices for a ring lattice.
        s = repelem((1:N)', 1, K);
        t = s + repmat(1:K, N, 1);
        t = mod(t-1, N) + 1;

        % Rewire the target node of each edge with probability beta
        for source = 1:N
            switchEdge = rand(K, 1) < beta;

            newTargets = rand(N, 1);
            newTargets(source) = 0;
            newTargets(s(t == source)) = 0;
            newTargets(t(source, ~switchEdge)) = 0;

            [~, ind] = sort(newTargets, 'descend');
            t(source, switchEdge) = ind(1:nnz(switchEdge));
        end

        h = graph(s, t);
end

function [width, height] = measureText(txt, textOpts, axis)
    if nargin < 3
       axis = gca(); 
    end
    if nargin < 2
        textOpts = struct();
        textOpts.HorizontalAlignment = 'center';
        textOpts.VerticalAlignment = 'middle';
        textOpts.FontSize = 20;
        textOpts.FontWeight = 'normal';
    end
    hTest = text(axis, 0, 0, txt, textOpts);
    textExt = get(hTest, 'Extent');
    delete(hTest);
    height = textExt(4)/3;    %Height
    width = textExt(3)/3;     %Width
end

function customeMarker_2(x, y, f)
    XY = [x, y];
    d = XY(1:end-1,:) - XY(2:end,:);
    slopex = d(:,2)./d(:,1);
    theta = atan(slopex); % Angle in radians
    theta = rad2deg(theta); % Convert to degrees
    p = XY(1:end-1,:) - 0.5*d;
    x0=p(:,1); y0=p(:,2);

    %ax = axes(f);
    % Triangle marker definition (relative size)
    markerSize = 0.08; % Keep this small so it doesn't scale with the figure
    %X_ = markerSize * [-1, 1, 0];  % X-coordinates (relative to center)
    %Y_ = markerSize * [-1, -1, 1]; % Y-coordinates 

    baseLength = 1;  % Length of the base
    height = 2;      % Height of the triangle
    X_ = markerSize * [-baseLength/2, baseLength/2, 0]; % X-coordinates (before rotation)
    Y_ = markerSize * [0, 0, height]; % Y-coordinates (before rotation)

    px = cell(length(x0),1);
    for k = 1:length(x0)
        if isnan(x0(k)), continue; end
        t = theta(k)-90;
        R = [cosd(t), -sind(t); sind(t), cosd(t)];
        rotatedXY = R * [X_; Y_]; % Apply rotation
        % Adjust for data point position
        %X = X_ + x0(k);
        %Y = Y_ + y0(k);        
        X = rotatedXY(1, :) + x0(k);
        Y = rotatedXY(2, :) + y0(k);        
        
        p = patch(X, Y, 'k', 'EdgeColor', 'k', 'LineWidth', .1);
        px{k} = p;

        %plot(x0, y0, 'k.', 'MarkerSize', 1); % Invisible anchor for reference        
        % addlistener(ax, 'MarkedClean', @(~,~) updatePatchSize(ax, p, x0(k), y0(k), markerSize));
    end
    % set(f, 'SizeChangedFcn', @(src,event) updatePatchSize(px, x0, y0, theta, markerSize));
    set(f, 'SizeChangedFcn', @(src,event) updatePatchSize_2(px, x, y, markerSize));
end




function customeMarker_1(x0, y0, theta, f)
    %ax = axes(f);
    % Triangle marker definition (relative size)
    markerSize = 0.08; % Keep this small so it doesn't scale with the figure
    %X_ = markerSize * [-1, 1, 0];  % X-coordinates (relative to center)
    %Y_ = markerSize * [-1, -1, 1]; % Y-coordinates 

    baseLength = 1;  % Length of the base
    height = 2;      % Height of the triangle
    X_ = markerSize * [-baseLength/2, baseLength/2, 0]; % X-coordinates (before rotation)
    Y_ = markerSize * [0, 0, height]; % Y-coordinates (before rotation)

    px = cell(length(x0),1);
    for k = 1:length(x0)
        if isnan(x0(k)), continue; end
        t = theta(k)-90;
        R = [cosd(t), -sind(t); sind(t), cosd(t)];
        rotatedXY = R * [X_; Y_]; % Apply rotation
        % Adjust for data point position
        %X = X_ + x0(k);
        %Y = Y_ + y0(k);        
        X = rotatedXY(1, :) + x0(k);
        Y = rotatedXY(2, :) + y0(k);        
        
        p = patch(X, Y, 'r', 'EdgeColor', 'k', 'LineWidth', .5);
        px{k} = p;

        %plot(x0, y0, 'k.', 'MarkerSize', 1); % Invisible anchor for reference        
        % addlistener(ax, 'MarkedClean', @(~,~) updatePatchSize(ax, p, x0(k), y0(k), markerSize));
    end
    set(f, 'SizeChangedFcn', @(src,event) updatePatchSize(px, x0, y0, theta, markerSize));
end


function customeMarker(x0, y0, theta)
    
    % Define triangle dimensions
    baseLength = .3;  % Length of the base
    height = .4;      % Height of the triangle
    % Define vertices relative to the base center
    X = [-baseLength/2, baseLength/2, 0]; % X-coordinates (before rotation)
    Y = [0, 0, height]; % Y-coordinates (before rotation)
    
    for k = 1:length(x0)
        if isnan(x0(k)), continue; end
        % Rotation matrix
        t = theta(k);        
        R = [cosd(t), -sind(t); sind(t), cosd(t)];
        rotatedXY = R * [X; Y]; % Apply rotation
        % Shift the rotated triangle to (x0, y0)
        X_rot = rotatedXY(1, :) + x0;
        Y_rot = rotatedXY(2, :) + y0;        
        patch(X_rot, Y_rot, 'r', 'EdgeColor', 'k', 'LineWidth', 1);
        %plot(x0, y0, 'k.', 'MarkerSize', 1); % Invisible anchor for reference
    end
end

% Callback function to update marker size when zooming
function updatePatchSize(patchObj, x0, y0, theta, markerSize)
    % fig = gcf;
    ax = gca;
    originalUnits = ax.Units; % Store original unit
    ax.Units = 'pixels';
    axPos = ax.Position;
    xLimits = xlim(ax);
    yLimits = ylim(ax);    
    dx = (xLimits(2) - xLimits(1)) * markerSize / axPos(3);
    dy = (yLimits(2) - yLimits(1)) * markerSize / axPos(4);
    ax.Units = originalUnits; % Restore original units (important!)

    dx = 1000*dx;
    dy = 1000*dy;

    baseLength = 1;  % Length of the base
    height = 2;      % Height of the triangle
    
    X_ = markerSize * [-baseLength/2, baseLength/2, 0]; % X-coordinates (before rotation)
    Y_ = markerSize * [0, 0, height]; % Y-coordinates (before rotation)
    
    fprintf('dx=%f dy=%f\n', dx, dy);

    for k = 1:length(x0)
        if isnan(x0(k)), continue; end
        t = theta(k)-90;
        R = [cosd(t), -sind(t); sind(t), cosd(t)];
        rotatedXY = R * [X_; Y_]; % Apply rotation

        X = rotatedXY(1, :) * dx + x0(k);
        Y = rotatedXY(2, :) * dy + y0(k);        
        
        %X = [-1 1 0] * dx + x0(k);
        %Y = [-1 -1 1] * dy + y0(k);
        set(patchObj{k}, 'XData', X, 'YData', Y);
    end
end


% Callback function to update marker size when zooming
function updatePatchSize_2(patchObj, x, y, markerSize)

    XY = [x, y];
    d = XY(1:end-1,:) - XY(2:end,:);
    slopex = d(:,2)./d(:,1);
    theta = atan(slopex); % Angle in radians
    theta = rad2deg(theta); % Convert to degrees
    p = XY(1:end-1,:) - 0.5*d;
    x0=p(:,1); y0=p(:,2);

    ax = gca;
    originalUnits = ax.Units; % Store original unit
    ax.Units = 'pixels';
    axPos = ax.Position;
    xLimits = xlim(ax);
    yLimits = ylim(ax);    
    dx = (xLimits(2) - xLimits(1)) * markerSize / axPos(3);
    dy = (yLimits(2) - yLimits(1)) * markerSize / axPos(4);
    ax.Units = originalUnits; % Restore original units (important!)

    dx = 1000*dx;
    dy = 1000*dy;

    baseLength = 1;       % Length of the base
    height = 2;           % Height of the triangle
    
    X_ = dx *  markerSize * [-baseLength/2, baseLength/2, 0]; % X-coordinates (before rotation)
    Y_ = dy * markerSize * [0, 0, height]; % Y-coordinates (before rotation)   
    
    for k = 1:length(x0)
        if isnan(x0(k)), continue; end
        t = theta(k)-90;
        R = [cosd(t), -sind(t); sind(t), cosd(t)];
        rotatedXY = R * [X_; Y_];
        X = rotatedXY(1, :) + x0(k);
        Y = rotatedXY(2, :) + y0(k);        
        set(patchObj{k}, 'XData', X, 'YData', Y);
    end
end
