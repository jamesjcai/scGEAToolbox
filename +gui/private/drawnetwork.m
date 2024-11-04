function [] = drawnetwork( x, y, nodeSizes, W, axis, varargin)
    isAxisSpecified = true;
    if(nargin < 5)
       isAxisSpecified = false;
    end
    p = inputParser;
    p.CaseSensitive = false;
%     addParameter(p, 'LineWidth', 0.75, @isnumeric);
    addParameter(p, 'XRatio', 1, @isnumeric);
    addParameter(p, 'OuterGap', 0.05, @isnumeric);
    addParameter(p, 'FontSize', 10, @isnumeric);
    addParameter(p, 'NodeLabels', {}, @iscell);
    addParameter(p, 'NodeColors', [], @isnumeric);
    addParameter(p, 'NodeCurvature', [], @isnumeric);
    addParameter(p, 'NodeFontSize', [], @isnumeric);
    addParameter(p, 'NodeLineWidth', [], @isnumeric);
    addParameter(p, 'NodeLineStyle', [], @iscell);
    addParameter(p, 'NodeLineColor', [], @isnumeric);
    addParameter(p, 'NodeLabelColor', [], @isnumeric);
    addParameter(p, 'EdgeColors', [], @isnumeric);
    addParameter(p, 'EdgeLineWidth', [], @isnumeric);
    addParameter(p, 'EdgeLineStyle', [], @iscell);
    addParameter(p, 'FitNodeSizes', [], @islogical);
    addParameter(p, 'DrawOuterRectangle', false, @islogical);
    parse(p, varargin{:});
    param = p.Results;

    nodeSizes = nodeSizes / 100;
    
    indices = find(triu(W, 1));
    [i1, i2] = ind2sub(size(W), indices);
    
%     linewidth = param.LineWidth;
    
    nNode = size(W, 1);
    nEdge = length(indices);
    
    if(~isempty(param.NodeFontSize))
        nodeFontSize = param.NodeFontSize;
    else
        nodeFontSize = repmat(param.FontSize, nNode, 1);
    end
    
    if(~isempty(param.NodeLabels))
        nodeLabels = param.NodeLabels;
    else
        nodeLabels = repmat({}, nNode, 1);
    end
    
%     nodeWidth = 2*nodeSizes;
    nodeWidth = nodeSizes;
    nodeHeight = nodeSizes;
%     nodeHeight = 2*nodeSizes;
    
    textOpts = struct();
    textOpts.HorizontalAlignment = 'center';
    textOpts.VerticalAlignment = 'middle';
    textOpts.FontSize = param.FontSize;
    textOpts.FontWeight = 'bold';
    
    if(~isempty(param.NodeLineWidth))
        nodeLineWidth = param.NodeLineWidth;
    else
        nodeLineWidth = repmat(0.5, nNode, 1);
    end
    
    gap = param.OuterGap;
    xratio = param.XRatio;
    aMin = min(x - nodeWidth/2)-gap;
    aMax = max(x + nodeWidth/2)+gap;
    bMin = min(y - nodeHeight/2)-gap;
    bMax = max(y + nodeHeight/2)+gap;
    aRange = aMax - aMin;
    bRange = bMax - bMin;
    cRange = max(aRange/xratio, bRange);

    xMin = aMin - (cRange*xratio - aRange)/2;
    xMax = aMax + (cRange*xratio - aRange)/2;
    yMin = bMin - (cRange - bRange)/2;
    yMax = bMax + (cRange - bRange)/2;
    
    xFactor = xMax - xMin;
    yFactor = yMax - yMin;
    
    if(~isAxisSpecified)
        axis = gca();
    end
    cla(axis);
    
    if(param.FitNodeSizes)
        hpadding = 0.05;
        vpadding = 0.01;
        xlim(axis, [0 1]);
        ylim(axis, [0 1]);
        for iNode = 1:nNode
            if(~isempty(nodeLabels{iNode}))
                linewidth = nodeLineWidth(iNode)/2;
                textOpts.FontSize = nodeFontSize(iNode);
                [width, height] = measureText(nodeLabels{iNode}, textOpts, axis);
                width = (width*(1+hpadding) + linewidth/100) * xFactor;
                height = (height*(1+vpadding) + linewidth/100) * yFactor;
                nodeWidth(iNode) = max(width, nodeWidth(iNode));
                nodeHeight(iNode) = max(height, nodeHeight(iNode));
            end
        end
    end
    
    if(~isempty(param.NodeCurvature))
        nodeCurvature = param.NodeCurvature;
    else
        nodeCurvature = ones(nNode, 1);
    end
    
    if(~isempty(param.NodeColors))
        nodeColors = param.NodeColors;
    else
        nodeColors = repmat([1 0.84 0], nNode, 1);
    end
    
    if(~isempty(param.NodeLineColor))
        nodeLineColor = param.NodeLineColor;
    else
        nodeLineColor = repmat([0 0 0], nNode, 1);
    end

    if(~isempty(param.NodeLineStyle))
        nodeLineStyle = param.NodeLineStyle;
    else
        nodeLineStyle = repmat({'-'}, nNode, 1);
    end
   
    if(~isempty(param.NodeLabelColor))
        nodeLabelColor = param.NodeLabelColor;
    else
        nodeLabelColor = repmat([0 0 0], nNode, 1);
    end
    
    if(~isempty(param.EdgeColors))
        edgeColors = param.EdgeColors;
    else
        edgeColors = repmat(0.65 * [1 1 1], nEdge, 1);
    end
    
    if(~isempty(param.EdgeLineWidth))
        edgeLineWidth = param.EdgeLineWidth;
    else
        edgeLineWidth = repmat(0.65 * [1 1 1], nEdge, 1);
    end

    if(~isempty(param.EdgeLineStyle))
        edgeLineStyle = param.EdgeLineStyle;
    else
        edgeLineStyle = repmat({'-'}, nEdge, 1);
    end
    
    rWidth = xMax - xMin;
    rHeight = yMax - yMin;
    gap = 0.002;
    rXMin = xMin + rWidth * gap;
    rXMax = xMax - rWidth * gap;
    rYMin = yMin + rHeight * gap;
    rYmax = yMax - rHeight * gap;
    
    hold(axis, 'on');
    rectangle(axis, 'Position', [rXMin rYMin (rXMax-rXMin) (rYmax-rYMin)], 'FaceColor', [1 1 1], 'EdgeColor', [0.7 0.7 0.7]);
%     rectangle(axis, 'Position', [-0.2 -0.2 1.4 1.4], 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
%     rectangle(axis, 'Position', [xMin yMin (xMax-xMin) (yMax-yMin)], 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
%     rectangle(axis, 'Position', [-0.1 -0.1 1.2 1.2], 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
    
%     tic
    groupEdges = true;
    if(groupEdges)
        [~, ~, ix] = unique(edgeLineStyle);
        [u, ib, ic] = unique([edgeLineWidth, edgeColors, ix], 'rows');
        nGroup = size(u, 1);
        if(nEdge == 1)
            ib = 1;
            ic = 1;
            nGroup= 1;
        end
        
        if(nEdge == 1)
           nGroup = 1; 
        end
        for iGroup = 1:nGroup
%         for iGroup = size(u, 1):-1:1
            ind = ic == iGroup;
            indices1 = i1(ind);
            indices2 = i2(ind);
            linewidth = edgeLineWidth(ib(iGroup));
            color = edgeColors(ib(iGroup), :);
            line_style = edgeLineStyle{ib(iGroup)};
            
            plot(axis, [x(indices1), x(indices2)]', ...
                [y(indices1), y(indices2)]', line_style, ...
                'Color', color, 'LineWidth', linewidth);
        end
%         plot(axis, [x(i1); x(i2)]', [y(i1); y(i2)]', '-', 'Color', [0.65 0.65 0.65]);
    else
        for i = 1:length(indices)
            linewidth = edgeLineWidth(i);
            ind1 = i1(i);
            ind2 = i2(i);
            x1 = x(ind1);
            y1 = y(ind1);
            x2 = x(ind2);
            y2 = y(ind2);
            color = edgeColors(i, :);
            plot(axis, [x1 x2], [y1 y2], '-', 'Color', color, 'LineWidth', linewidth);
        end
    end
%     toc
    
    for i = 1:nNode
        width = nodeWidth(i);
        height = nodeHeight(i);
        xx = x(i) - width/2;
        yy = y(i) - height/2;
%         curvature = 1;
        color = nodeColors(i, :);
        curvature = nodeCurvature(i);
        nodelinewidth = nodeLineWidth(i);
        nodelinecolor = nodeLineColor(i, :);
        nodelabelcolor = nodeLabelColor(i, :);
        nodelinestyle = nodeLineStyle{i};
        rectangle(axis, 'Position', [xx yy width height], ...
            'Curvature', curvature, ...
            'FaceColor', color, ...
            'LineWidth', nodelinewidth, ...
            'LineStyle', nodelinestyle, ...
            'EdgeColor', nodelinecolor);
        if(~isempty(nodeLabels{i}))
            txt = nodeLabels{i};
            textOpts.FontSize = nodeFontSize(i);
            text(axis, x(i), y(i), txt, 'Color', nodelabelcolor, textOpts); 
        end
    end
    hold(axis, 'off');
    
    if(param.DrawOuterRectangle)
        rWidth = xMax - xMin;
        rHeight = yMax - yMin;
        gap = 0.005;
        rXMin = xMin + rWidth * gap;
        rXMax = xMax - rWidth * gap;
        rYMin = yMin + rHeight * gap;
        rYmax = yMax - rHeight * gap;
        
%         rectangle(axis, ...
%             'Position', [rXMin rYMin (rXMax-rXMin) (rYmax-rYMin)], ...
%             'FaceColor', 'none', ...
%             'EdgeColor', [0.7 0.7 0.7], ...
%             'LineWidth', 0.75);
    end
    
    xlim(axis, [xMin xMax]);
    ylim(axis, [yMin yMax]);
    set(axis, 'Visible', 'off');
end

function[width, height] = measureText(txt, opt, axis)
    if(nargin < 3)
       axis = gca(); 
    end
    hTest = text(axis, 0, 0, txt, opt);
    textExt = get(hTest, 'Extent');
    delete(hTest);
    height = textExt(4)/3;    %Height
    width = textExt(3)/3;     %Width
end
