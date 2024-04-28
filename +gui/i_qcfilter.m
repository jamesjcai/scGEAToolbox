function [X, idx] = i_qcfilter(X, genelist)
% https://www.mathworks.com/matlabcentral/answers/143306-how-to-move-a-plotted-line-vertically-with-mouse-in-a-gui

nGenes = sum(X > 0, 1)';
%nUMIs=sum(X)';
mtratio = sc_mtratio(X, genelist);
xr = [300, 5000];
yr = [0.025, 0.1];

x = nGenes;
y = mtratio;

fh = figure('DockControls', 'off');
scatter(x, y);
xlim([min([xr(1) * 0.9, min(x) * 0.9]), max([xr(2) * 1.1, max(x) * 1.1])]);
ylim([min([yr(1) * 0.9, min(y) * 0.9]), max([yr(2) * 1.1, max(y) * 1.1])]);
lh1 = xline(xr(1), 'r-');
lh2 = xline(xr(2), 'r-');
lh3 = yline(yr(1), 'r-');
lh4 = yline(yr(2), 'r-');
idx = true(length(nGenes), 1);
guidata(fh, [x, y]);
set(fh, 'WindowButtonDownFcn', @mouseDownCallback);
i = (x > lh1.Value) & (x < lh2.Value);
j = (y > lh3.Value) & (y < lh4.Value);
idx = i & j;
title(sprintf('%d out of %d (%.2f%%)', ...
    sum(idx), length(i), 100*sum(idx)./length(idx)));
uiwait(fh);
X = X(:, idx);

if ~(ismcc || isdeployed)
    labels = {'Save X to variable named:', ...
        'Save IDX to variable named:'};
    vars = {'X', 'idx'};
    values = {X, idx};
    export2wsdlg(labels, vars, values, ...
        'Save Data to Workspace', ...
        logical([1, 0]));
end

% figure;
% subplot(2,2,1)
% scatter(nGenes,mtratio);
% xlabel('nGenes'); ylabel('pMito');
% xline(300,'r-'); xline(5000,'r-');
% yline(0.025,'r-'); yline(0.1,'r-');
%
% subplot(2,2,2)
% scatter(nUMIs,mtratio);
% xlabel('nUMI'); ylabel('pMito');
% xline(400,'r-'); xline(30000,'r-');
% yline(0.025,'r-'); yline(0.1,'r-');
%
%
% subplot(2,2,3)
% scatter(nGenes,nUMIs);
% xlabel('nGenes'); ylabel('nUMI');
% xline(300,'r-'); xline(5000,'r-');
% yline(400,'r-'); yline(30000,'r-');

% https://github.com/prabhakarlab/RCAv2


    function mouseDownCallback(figHandle, varargin)
        xydata = guidata(figHandle);
        currentPoint = get(figHandle, 'CurrentPoint');
        x1 = currentPoint(1, 1);
        y1 = currentPoint(1, 2);
        axes1 = get(figHandle, 'CurrentAxes');
        set(axes1, 'Units', 'pixels');
        axesPos = get(axes1, 'Position');
        minx = axesPos(1);
        miny = axesPos(2);
        maxx = minx + axesPos(3);
        maxy = miny + axesPos(4);
        if x1 >= minx && x1 <= maxx && y1 >= miny && y1 <= maxy
            currentPoint = get(axes1, 'CurrentPoint');
            xx = currentPoint(2, 1);
            yy = currentPoint(2, 2);
            if min(abs(xx-axes1.XLim)/diff(axes1.XLim)) < min(abs(yy-axes1.YLim)/diff(axes1.YLim))
                if abs(xx-lh1.Value) < abs(xx-lh2.Value)
                    if ~isempty(lh1), delete(lh1); end
                    lh1 = xline(xx, 'r-');
                else
                    if ~isempty(lh2), delete(lh2); end
                    lh2 = xline(xx, 'r-');
                end
            else
                if abs(yy-lh3.Value) < abs(yy-lh4.Value)
                    if ~isempty(lh3), delete(lh3); end
                    lh3 = yline(yy, 'r-');
                else
                    if ~isempty(lh4), delete(lh4); end
                    lh4 = yline(yy, 'r-');
                end
            end
            i = (xydata(:, 1) > lh1.Value) & (xydata(:, 1) < lh2.Value);
            j = (xydata(:, 2) > lh3.Value) & (xydata(:, 2) < lh4.Value);
            idx = i & j;
            xr = [lh1.Value, lh2.Value];
            yr = [lh3.Value, lh4.Value];
            axes1.Title.String = sprintf('%d out of %d (%.2f%%)', ...
                sum(i & j), length(i), 100*sum(idx)./length(i));
        end
end


end
