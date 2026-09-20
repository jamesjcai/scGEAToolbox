function i_heatmap(sce, glist, thisc, parentfig)

if nargin<4, parentfig = []; end

[c, cL, noanswer] = gui.i_reordergroups(thisc, [], parentfig);
if noanswer, return; end
[~, gidx] = ismember(glist, sce.g);
[Xt] = gui.i_transformx(sce.X, true, 3, parentfig);
if isempty(Xt), return; end

Y = Xt(gidx, :);
[~, cidx] = sort(c);
Yori = Y(:, cidx);

[methodid, dim] = gui.i_selnormmethod(parentfig);
if isempty(dim) || isempty(methodid), return; end
[Y] = gui.i_norm4heatmap(Yori, dim, methodid);

% szgn = grpstats(c, c, @numel);
szgn = splitapply(@numel, c, c);
a = zeros(1, max(c));
b = zeros(1, max(c));
for kx = 1:max(c)
    a(kx) = sum(c <= kx);
    b(kx) = round(sum(c == kx)./2);
end

% figure;
% heatmap(Y)
% assignin('base','Y',Y);
% assignin('base','g',glist) ;
% heatmap(Y,'YDisplayLabels',glist, ...
%     'XDisplayLabels',strings(size(Y,2),1), ...
%     'GridVisible',false,'ColorScaling','scaled',...
%     'ColorbarVisible',false)

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

% Where each column block sits now, in terms of the groups GUI.I_REORDER-
% GROUPS handed over. IN_CALLBACK_SORTGROUPS permutes the display and
% keeps this in step, so "Unsorted" has an order to go back to.
ordnow = (1:length(cL)).';

% for k = 1:length(cL)
%     Z(:, k) = mean(Y(:, c == k), 2);
% end

% [Z] = gui.i_norm4heatmap(Z);

% figure;
% h2=heatmap(strrep(cL,'_','\_'),MX,Z);
% h2.Title = 'Marker Gene Heatmap';
% h2.XLabel = 'Group';
% h2.YLabel = 'Marker Gene';
% h2.Colormap = parula;
% h2.GridVisible = 'off';
% h2.CellLabelColor='none';
% h2.ColorLimits=[min(Z(:)), max(Z(:))];

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
        h = gui.i_drawgroupmap(gca, h, drawspec);
    end

function in_callback_renamecat(~, ~)
        tg = gui.i_inputgenelist(string(cL), true, hFig);
        if isempty(tg), return; end
        if length(tg) == length(cL)
            set(gca, 'XTick', a-b);
            set(gca, 'XTickLabel', tg(:))
            cL = tg;
        else
            gui.myErrordlg(hFig, 'Wrong input.');
        end
    end

function in_callback_resetcolor(~, ~)
        set(gca, 'FontSize', 10);
        colormap default
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

function in_callback_exporttable(~, ~, T, needwait, defname)
        if nargin < 5, defname = []; end
        if nargin < 4, needwait = false; end
        if ~isempty(defname)
            [file, path] = uiputfile({'*.txt'; '*.*'}, 'Save as', defname);
        else
            [file, path] = uiputfile({'*.txt'; '*.*'}, 'Save as');
        end
        if pkg.i_isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure'), figure(parentfig); end
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
                gui.myHelpdlg(hFig, ...
                    sprintf('Result has been saved in %s', filename));
            end
        end
    end


function in_callback_summarymap(~, ~)
        for ky = 1:length(cL)
            Z(:, ky) = mean(Y(:, c == ky), 2);
        end

        hx1=gui.myFigure(parentfig);

        [mx,idx]=unique(MX,'stable');
        z = Z(idx,:);
        h = heatmap(gui.i_escapeunderscore(cL), mx, z);
        h.Title = 'Marker Gene Heatmap';
        h.XLabel = 'Group';
        h.YLabel = 'Marker Gene';
        h.Colormap = parula;
        h.GridVisible = 'off';
        h.CellLabelColor = 'none';
        t = array2table(z, 'VariableNames', cL, 'RowNames', mx);
        % writetable(t,'aaa.csv','WriteRowNames',true);
        hx1.addCustomButton('off', {@in_callback_exporttable, t}, 'floppy-disk-arrow-in.jpg', 'Save table...');
        hx1.addCustomButton('off', @in_callback_resetcolor, 'refresh_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Reset color map');
        % disp('https://software.broadinstitute.org/morpheus/')
        hx1.show(hFig);
    end

function in_callback_summarymapT(~, ~)
            for ky = 1:length(cL)
                Z(:, ky) = mean(Y(:, c == ky), 2);
            end

        hx2=gui.myFigure(parentfig);
        % assignin("base","Z",Z);
        % assignin("base","MX",MX);
        % assignin("base","cL",cL);
        % matlab.lang.makeUniqueStrings(MX)

        [mx,idx]=unique(MX,'stable');
        z = Z(idx,:);
        h = heatmap(mx, gui.i_escapeunderscore(cL), z.');
        h.Title = 'Marker Gene Heatmap';
        h.YLabel = 'Group';
        h.XLabel = 'Marker Gene';
        h.Colormap = parula;
        h.GridVisible = 'off';
        h.CellLabelColor = 'none';
        t = array2table(z.', 'VariableNames', mx, 'RowNames', cL);
        %         s = struct(h);
        %         s.XAxis.TickLabelRotation=45;
        % writetable(t,'aaa.csv','WriteRowNames',true);
        hx2.addCustomButton('off', {@in_callback_exporttable, t}, 'floppy-disk-arrow-in.jpg', 'Save table...');
        hx2.addCustomButton('off', {@gui.i_pickcolormap, c}, 'color-wheel.jpg', 'Pick new color map...');
        hx2.addCustomButton('off', @in_callback_resetcolor, 'refresh_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Reset color map');
        hx2.show(hFig);
    end

function in_callback_dotplotx(~, ~)
        try
            gui.i_dotplot(sce.X, sce.g, c, cL, MX);
        catch ME
            gui.myErrordlg(hFig, ME.message, ME.identifier);
        end
    end

end
