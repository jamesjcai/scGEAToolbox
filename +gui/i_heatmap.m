function i_heatmap(sce, glist, thisc, parentfig)

if nargin<4, parentfig = []; end

[c, cL, noanswer] = gui.i_reordergroups(thisc, [], parentfig);
if noanswer, return; end
[~, gidx] = ismember(glist, sce.g);
[Xt] = gui.i_transformx(sce.X, true, "libsize_log1p", parentfig);
if isempty(Xt), return; end

Y = Xt(gidx, :);
[~, cidx] = sort(c);
Yori = Y(:, cidx);

[methodid, dim] = gui.i_selnormmethod(parentfig);
if isempty(dim) || isempty(methodid), return; end
[Y] = gui.i_norm4heatmap(Yori, dim, methodid);

szgn = splitapply(@numel, c, c);
a = zeros(1, max(c));
b = zeros(1, max(c));
for kx = 1:max(c)
    a(kx) = sum(c <= kx);
    b(kx) = round(sum(c == kx)./2);
end

hx=gui.myFigure(parentfig);
hFig=hx.FigHandle;

% The first draw goes through IN_DRAWMAP like every later one, so the map
% the user is handed and the map a button leaves behind cannot drift apart.
% That is how the group separators came to be lost: the lines were drawn
% once, here, and every redraw below emptied the axes without knowing to
% put them back. H has to exist for the DELETE at the top of IN_DRAWMAP;
% an empty handle array deletes to nothing.
h = gobjects(0);
fliped = false;
in_drawmap();

hx.addCustomButton('off', @in_callback_renamecat, 'edit.jpg', 'Rename groups...');
hx.addCustomButton('off', @in_callback_sortgroups, 'reorder.jpg', 'Sort groups...');
hx.addCustomButton('off', @in_callback_resetcolor, 'refresh_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Reset color map');
hx.addCustomButton('off', @in_callback_flipxy, 'mat-wrap-text.jpg', 'Flip XY');
hx.addCustomButton('on', @in_callback_summarymap, 'HDF_object01.gif', 'Summary map...');
hx.addCustomButton('off', @in_callback_summarymapT, 'HDF_object02.gif', 'Summary map, transposed...');
hx.addCustomButton('on', @in_callback_savetable, 'floppy-disk-arrow-in.jpg', 'Export data...');
hx.addCustomButton('off', @in_callback_changenorm, 'mw-pickaxe-mining.jpg', 'Change normalization method...');
hx.addCustomButton('off', @in_callback_dotplotx, 'icon-mat-blur-linear-10.gif', 'Dot plot...');

hx.show(parentfig);

MX = glist;

c = c(cidx);
Z = zeros(length(glist), length(cL));

% The columns on screen, as indices into SCE's cells. C is in drawn order,
% so anything that goes back to SCE.X has to take its cells in this order
% too; IN_CALLBACK_SORTGROUPS permutes it with YORI.
colorder = cidx(:).';

% Where each column block sits now, in terms of the groups GUI.I_REORDER-
% GROUPS handed over. IN_CALLBACK_SORTGROUPS permutes the display and
% keeps this in step, so "Unsorted" has an order to go back to.
ordnow = (1:length(cL)).';

function in_callback_changenorm(~, ~)
        % HFIG, not PARENTFIG: the button is on the heatmap figure, and
        % every other callback here parents its dialogs to HFIG.
        % GUI.I_SELNORMMETHOD calls FIGURE() on what it is given, so
        % PARENTFIG raised the main scgeatool window over the heatmap the
        % user had just clicked in, and centred the dialog there.
        [methodid, dim] = gui.i_selnormmethod(hFig);

        % Cancelling either dialog leaves both empty, and the SWITCH
        % METHODID in GUI.I_NORM4HEATMAP threw on the empty rather than
        % the heatmap being left as it was. The call at the top of the
        % file already returns on this; this one redrew.
        if isempty(dim) || isempty(methodid), return; end

        [Y] = gui.i_norm4heatmap(Yori, dim, methodid);
        % Y = log1p(Y);
        % Through IN_DRAWMAP rather than drawing here. This drew the map
        % the upright way round whatever FLIPED said, so renormalizing a
        % flipped map quietly un-flipped it while the button still thought
        % it was flipped, and it lost the group separators the same way
        % the flip did.
        in_drawmap();
    end

function in_callback_flipxy(~, ~)
        % This used to transpose the image and move the ticks, and draw no
        % group separators. IMAGESC empties the axes it draws into, so the
        % yellow lines marking where one group ends and the next begins
        % went with the old image and never came back: one click and the
        % map had no boundaries on it for the rest of its life, in either
        % orientation. IN_DRAWMAP draws the map whole, lines included.
        fliped = ~fliped;
        in_drawmap();
    end

function in_callback_sortgroups(~, ~)
        % Reorder the column blocks. Only the order changes: the values,
        % the genes and the normalization all stay where they are, so this
        % never has to recompute anything.
        %
        % The permutation is worked out against the order on screen and
        % applied to it, rather than rebuilt from the original each time,
        % which is what lets a sort compose with a rename.
        sortby = gui.i_askgrouporder(hFig);
        if sortby == "", return; end

        [newc, p, neworder] = gui.i_grouporderperm(c, cL, ordnow, sortby);
        if isequal(p, (1:numel(cL)).'), return; end

        c = newc;
        cL = cL(p);
        ordnow = ordnow(p);
        Yori = Yori(:, neworder);
        Y = Y(:, neworder);
        colorder = colorder(neworder);

        szgn = splitapply(@numel, c, c);
        a = zeros(1, max(c));
        b = zeros(1, max(c));
        for kg = 1:max(c)
            a(kg) = sum(c <= kg);
            b(kg) = round(sum(c == kg)./2);
        end
        in_drawmap();
    end

function in_drawmap()
        % The one place the map is drawn, so the first draw and every
        % redraw agree. GUI.I_DRAWGROUPMAP does the work and says why it
        % is one function rather than a branch in each callback.
        drawspec = struct( ...
            'Y',      Y, ...
            'cL',     {cL}, ...
            'glist',  {glist}, ...
            'a',      a, ...
            'b',      b, ...
            'szgn',   szgn, ...
            'fliped', fliped);
        h = gui.i_drawgroupmap(hx.AxHandle, h, drawspec);
    end

function in_callback_renamecat(~, ~)
        tg = gui.i_inputgenelist(string(cL), true, hFig);
        if isempty(tg), return; end
        if length(tg) == length(cL)
            % Redrawn rather than relabelled in place: the groups are on
            % the Y axis once the map is flipped, and setting XTickLabel
            % wrote the new names over the gene labels there.
            cL = tg;
            in_drawmap();
        else
            gui.myErrordlg(hFig, 'Wrong input.');
        end
    end

function in_callback_resetcolor(src, ~)
        % Shared by the main map and the two summary maps, so the target is
        % the clicked window's plot, not gca -- which on a summary window
        % would lay a new empty axes over the chart. A heatmap chart keeps
        % its own Colormap and starts from parula.
        fig = ancestor(src, 'figure');
        target = fig.CurrentAxes;
        if isempty(target), return; end
        set(target, 'FontSize', 10);
        if isa(target, 'matlab.graphics.chart.HeatmapChart')
            target.Colormap = parula;
        else
            colormap(target, 'default');
        end
    end


function in_callback_savetable(~, ~)
        labels = {'Save Y to variable named:', ...
            'Save glist to variable named:', ...
            'Save cL to variable named:'};
        vars = {'Y', 'g', 'cL'};
        values = {full(Y), glist, string(cL)};
        [~, ~] = export2wsdlg(labels, vars, values, ...
            'Save Data to Workspace');
    end

function in_callback_exporttable(src, ~, T, needwait, defname)
        % Used by the summary maps: after the file dialog, raise the window
        % whose Save button was clicked, not the main app (PARENTFIG).
        srcfig = ancestor(src, 'figure');
        if nargin < 5, defname = []; end
        if nargin < 4, needwait = false; end
        if ~isempty(defname)
            [file, path] = uiputfile({'*.txt'; '*.*'}, 'Save as', defname);
        else
            [file, path] = uiputfile({'*.txt'; '*.*'}, 'Save as');
        end
        gui.i_raisefig(srcfig);
        if isequal(file, 0) || isequal(path, 0)
            return;
        else
            filename = fullfile(path, file);
            try
                writetable(T, filename, 'Delimiter', '\t', 'WriteRowNames', true);
            catch
                writematrix(T, filename, 'Delimiter', '\t');
            end
            drawnow;
            if needwait
                gui.myHelpdlg(srcfig, ...
                    sprintf('Result has been saved in %s', filename));
            end
        end
    end


function in_callback_summarymap(~, ~)
        gui.myFigure.drawInto(hFig, @in_summarymap);
    end

function in_summarymap()
        for ky = 1:length(cL)
            Z(:, ky) = mean(Y(:, c == ky), 2);
        end

        hx1=gui.myFigure(parentfig);

        [mx,idx]=unique(MX,'stable');
        z = Z(idx,:);
        % HS, not H: these callbacks share the parent's workspace, and H
        % is the main map's image, which IN_DRAWMAP deletes on the next
        % flip, sort or renormalization -- taking this chart with it.
        hs = heatmap(hx1.FigHandle, gui.i_escapeunderscore(cL), mx, z);
        hs.Title = 'Marker Gene Heatmap';
        hs.XLabel = 'Group';
        hs.YLabel = 'Marker Gene';
        hs.Colormap = parula;
        hs.GridVisible = 'off';
        hs.CellLabelColor = 'none';
        t = array2table(z, 'VariableNames', cL, 'RowNames', mx);
        hx1.addCustomButton('off', {@in_callback_exporttable, t}, 'floppy-disk-arrow-in.jpg', 'Save table...');
        hx1.addCustomButton('off', @in_callback_resetcolor, 'refresh_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Reset color map');
        hx1.show(hFig);
    end

function in_callback_summarymapT(~, ~)
        gui.myFigure.drawInto(hFig, @in_summarymapT);
    end

function in_summarymapT()
            for ky = 1:length(cL)
                Z(:, ky) = mean(Y(:, c == ky), 2);
            end

        hx2=gui.myFigure(parentfig);

        [mx,idx]=unique(MX,'stable');
        z = Z(idx,:);
        hs = heatmap(hx2.FigHandle, mx, gui.i_escapeunderscore(cL), z.');
        hs.Title = 'Marker Gene Heatmap';
        hs.YLabel = 'Group';
        hs.XLabel = 'Marker Gene';
        hs.Colormap = parula;
        hs.GridVisible = 'off';
        hs.CellLabelColor = 'none';
        t = array2table(z.', 'VariableNames', mx, 'RowNames', cL);
        hx2.addCustomButton('off', {@in_callback_exporttable, t}, 'floppy-disk-arrow-in.jpg', 'Save table...');
        hx2.addCustomButton('off', {@gui.i_pickcolormap, c}, 'color-wheel.jpg', 'Pick new color map...');
        hx2.addCustomButton('off', @in_callback_resetcolor, 'refresh_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Reset color map');
        hx2.show(hFig);
    end

function in_callback_dotplotx(~, ~)
        try
            % C is in drawn order, so the cells must be too: SCE.X in its
            % own order gave every group another group's cells.
            gui.myFigure.drawInto(hFig, ...
                @() gui.i_dotplot(sce.X(:, colorder), sce.g, c, cL, MX));
        catch ME
            gui.myErrordlg(hFig, ME.message, ME.identifier);
        end
    end

end
