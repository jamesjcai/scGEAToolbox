function [h1, h2] = sc_scattermarker(X, genelist, ...
    s, targetg, methodid, ...
    sz, showcam)
% SC_SCATTERMARKER(X,genelist,g,s,methodid)
%
% USAGE:
% s=sc_tsne(X,3);
% g=["AGER","SFTPC","SCGB3A2","TPPP3"];
% sc_scattermarker(X,genelist,s,genelist(1));
%
% INPUTS:
% X: Numeric matrix of gene expression data (genes x cells).
% genelist: Cell array or string array of gene names (length = rows of X).
% s: Matrix of 2D or 3D coordinates for cells (cells x 2 or cells x 3).
% targetg: Gene(s) to visualize (string, char, or cell array).
% methodid: (Optional) Visualization method (1 to 5).
% sz: (Optional) Marker size for scatter plots (default = 5).
% showcam: (Optional) Enable/disable camera controls (default = true).
%
% OUTPUTS:
% h1: Handle to the primary axes.
% h2: Handle to the secondary axes (if applicable).

import gui.*
h1 = [];
h2 = [];
if nargin < 4
    error('sc_scattermarker(X, genelist, s, g)');
end

if ~isnumeric(X) || ~ismatrix(X)
    error('Input X must be a numeric matrix.');
end

if size(X, 1) ~= length(genelist)
    error('The number of rows in X must match the length of genelist.');
end


if isvector(s) || isscalar(s)
    error('S should be a matrix.');
end
if nargin < 7, showcam = true; end
if nargin < 6 || isempty(sz), sz = 5; end
if nargin < 5, methodid = 1; end

% Normalize targetg to string array for uniform handling
if ischar(targetg)
    targetg = string(targetg);
elseif iscellstr(targetg)
    targetg = string(targetg);
end

if ~isStringScalar(targetg)
    for k = 1:length(targetg)
        figure;
        sc_scattermarker(X, genelist, s, targetg(k), methodid, sz);
    end
elseif isStringScalar(targetg)
    if ismember(targetg, genelist)
        x = s(:, 1);
        y = s(:, 2);
        if min(size(s)) == 2
            z = [];
        else
            z = s(:, 3);
        end
        idx = strcmp(genelist, targetg);
        c = full(X(idx, :));

        titxt = '';
        switch methodid
            case 1
                gui.i_stemscatter([x, y], c);
                h1 = gca;
                grid(h1, "on");
                title(h1, targetg);
                titxt = gui.i_getsubtitle(c);
                subtitle(h1, titxt);
            case 2
                if isempty(z)
                    scatter(x, y, sz, c, 'filled');
                else
                    scatter3(x, y, z, sz, c, 'filled');
                end
                h1 = gca;
                set(h1, 'XTickLabel', []);
                set(h1, 'YTickLabel', []);
                set(h1, 'ZTickLabel', []);
                grid(h1, "on");

                title(h1, targetg);

                titxt = gui.i_getsubtitle(c);
                subtitle(h1, titxt);

            case 3

            case 4
                h1 = subplot(1, 2, 1);
                sc_scattermarker(X, genelist, s, targetg, 2, sz, false);
                h2 = subplot(1, 2, 2);
                sc_scattermarker(X, genelist, s, targetg, 1, sz, false);
                hFig = gcf;
                hFig.Position(3) = hFig.Position(3) * 2;
            case 5 % ============ 5

                if size(s, 2) >= 3
                    x = s(:, 1);
                    y = s(:, 2);
                    z = s(:, 3);
                    is2d = false;
                else
                    x = s(:, 1);
                    y = s(:, 2);
                    z = zeros(size(x));
                    is2d = true;
                end

                h1 = subplot(1, 2, 1);
                s1 = scatter3(h1, x, y, z, sz, c, 'filled');
                set(h1, 'XTickLabel', []);
                set(h1, 'YTickLabel', []);
                set(h1, 'ZTickLabel', []);
                grid(h1, 'on');
                title(h1, targetg);

                titxt = gui.i_getsubtitle(c);
                subtitle(h1, titxt);

                h2 = subplot(1, 2, 2);
                s2 = stem3(h2, x, y, c, 'marker', 'none', 'color', 'm');
                hold(h2, 'on');
                scatter3(h2, x, y, zeros(size(y)), 5, c, 'filled');
                %[ayy,byy]=view(h2);
                set(h2, 'XTickLabel', []);
                set(h2, 'YTickLabel', []);
                grid(h2, 'on');

                if ~is2d
                    evalin('base', 'linkprop(findobj(gcf,''type'',''axes''), {''CameraPosition'',''CameraUpVector''});');
                    evalin('base', 'rotate3d on');
                end
                hFig = gcf;
                hFig.Position(3) = hFig.Position(3) * 1.8;
                if ~is2d
                    view(h1, 3);
                else
                    view(h2, 3);
                end
        otherwise
                warning('sc_scattermarker:InvalidMethod', 'Unknown methodid %d. Use 1–5.', methodid);
        end % end of switch of methodid


        a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');
        gui.i_setautumncolor(c, a, true, any(c == 0), [], hFig);

        ori_c = c;

        if showcam
            hFig = gcf;
            tb = uitoolbar(hFig);
            gui.i_addbutton2fig(tb, 'off', @gui.i_linksubplots, 'keyframes-minus.jpg', 'Link subplots');
            gui.i_addbutton2fig(tb, 'on', {@i_genecards, targetg}, 'www.jpg', 'GeneCards...');
            gui.i_addbutton2fig(tb, 'on', {@i_PickColorMap, c}, 'plotpicker-compass.gif', 'Pick new color map...');
            gui.i_addbutton2fig(tb, 'off', @i_RescaleExpr, 'IMG00074.GIF', 'Rescale expression level [log2(x+1)]');
            gui.i_addbutton2fig(tb, 'off', @i_ResetExpr, 'plotpicker-geobubble2.gif', 'Reset expression level');
            gui.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
        end
    else
        warning('%s no expression', targetg);
    end

    if showcam
        gui.i_3dcamera(tb, targetg, false, hFig);
    end
end

    function i_RescaleExpr(~, ~)
    c = log2(1+c);
    [ax, bx] = view(h2);
    delete(s2);
    s2 = stem3(h2, x, y, c, 'marker', 'none', 'color', 'm');
    view(h2, ax, bx);
    title(h2, targetg);
    subtitle(h2, titxt);

    [ax, bx] = view(h1);
    delete(s1);
    s1 = scatter3(h1, x, y, z, sz, c, 'filled');
    view(h1, ax, bx);
    colorbar(h1);

    title(h1, targetg);
    subtitle(h1, titxt);
    end

    function i_ResetExpr(~, ~)
    c = ori_c;
    delete(s2);
    s2 = stem3(h2, x, y, c, 'marker', 'none', 'color', 'm');

    delete(s1);
    s1 = scatter3(h1, x, y, z, sz, c, 'filled');
    title(h1, targetg);
    subtitle(h1, titxt);
    title(h2, targetg);
    subtitle(h2, titxt);
    end

end

    function i_genecards(~, ~, g)
        web(sprintf('https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', g));
    end

    function i_PickColorMap(src, ~, c)
        % A local function: hFig is not in scope here, so use the button's figure.
        hFig = ancestor(src, 'figure');
        cleanupObj = onCleanup(@() gui.i_raisefig(hFig));
        list = {'parula', 'turbo', 'hsv', 'hot', 'cool', 'spring', ...
            'summer', 'autumn (default)', ...
            'winter', 'jet'};
        [indx, tf] = listdlg('ListString', list, 'SelectionMode', 'single', ...
            'PromptString', 'Select a colormap:', 'ListSize', [220, 300]);
        if tf == 1
            a = list{indx};
            if strcmp(a, 'autumn (default)')
                a = 'autumn';
            end
            gui.i_setautumncolor(c, a, [], [], hFig.CurrentAxes, hFig);
            setpref('scgeatoolbox', 'prefcolormapname', a);
        end
    end



    function selectcolormapeditor(~, ~)
        % colormapeditor;
    end

    function within_stemscatter(x, y, z)
        if nargin < 3
            x = randn(300, 1);
            y = randn(300, 1);
            z = abs(randn(300, 1));
        end
        if isempty(z)
            warning('No expression');
            scatter(x, y, '.');
        else
            stem3(x, y, z, 'marker', 'none', 'color', 'm');
            hold on;
            scatter(x, y, 10, z, 'filled');
            hold off;
        end
    end
