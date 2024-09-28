function callback_Heatmap(src, ~)
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
[thisc, ~] = gui.i_select1class(sce);
if isempty(thisc), return; end

[c, cL, noanswer] = gui.i_reordergroups(thisc);
if noanswer, return; end

[glist] = gui.i_selectngenes(sce,[],FigureHandle);
if isempty(glist)
    helpdlg('No gene selected.', '');
    return;
end

[~, gidx] = ismember(glist, sce.g);

%[Xt]=gui.i_transformx(sce.X);
%if isempty(Xt), return; end
Xt = sc_norm(sce.X);
Xt = log1p(Xt);

Y = Xt(gidx, :);
[~, cidx] = sort(c);
Yori = Y(:, cidx);
[Y] = gui.i_norm4heatmap(Yori);

szgn = grpstats(c, c, @numel);
a = zeros(1, max(c));
b = zeros(1, max(c));
for k = 1:max(c)
    a(k) = sum(c <= k);
    b(k) = round(sum(c == k)./2);
end

% figure;
% heatmap(Y)
% assignin('base','Y',Y);
% assignin('base','g',glist) ;
% heatmap(Y,'YDisplayLabels',glist, ...
%     'XDisplayLabels',strings(size(Y,2),1), ...
%     'GridVisible',false,'ColorScaling','scaled',...
%     'ColorbarVisible',false)

hFig = figure('Visible', 'off');
h = imagesc(Y);
% hFig.Colormap = repmat(linspace(0, 1, 25).', 1, 3);
set(gca, 'XTick', a-b);
set(gca, 'XTickLabel', strrep(cL, '_', '\_'));
%set(gca,'XTickLabelRotation',0);
set(gca, 'YTick', 1:length(glist));
set(gca, 'YTickLabel', glist);
set(gca, 'TickLength', [0, 0]);
% colormap(flipud(bone));
box on

szc = cumsum(szgn);
for k = 1:length(szc)
    xline(szc(k)+0.5, 'y-');
end

tb = findall(hFig, 'Tag', 'FigureToolBar');
uipushtool(tb, 'Separator', 'off');
pkg.i_addbutton2fig(tb, 'on', {@gui.i_pickcolormap, c}, 'plotpicker-compass.gif', 'Pick new color map...');
pkg.i_addbutton2fig(tb, 'off', @gui.i_changefontsize, 'noun_font_size_591141.gif', 'ChangeFontSize');
pkg.i_addbutton2fig(tb, 'on', @i_renamecat, 'guideicon.gif', 'Rename groups...');
pkg.i_addbutton2fig(tb, 'on', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
pkg.i_addbutton2fig(tb, 'on', @gui.i_invertcolor, 'plotpicker-comet.gif', 'Invert colors');
pkg.i_addbutton2fig(tb, 'off', @i_resetcolor, 'plotpicker-geobubble2.gif', 'Reset color map');
pkg.i_addbutton2fig(tb, 'off', @i_flipxy, 'xplotpicker-geobubble2.gif', 'Flip XY');
pkg.i_addbutton2fig(tb, 'off', @i_summarymap, 'HDF_object01.gif', 'Summary map...');
pkg.i_addbutton2fig(tb, 'off', @i_summarymapT, 'HDF_object02.gif', 'Summary map, transposed...');
pkg.i_addbutton2fig(tb, 'on', {@gui.i_resizewin, hFig}, 'HDF_pointx.gif', 'Resize Plot Window');


try
    [px_new] = gui.i_getchildpos(FigureHandle, hFig);
    if ~isempty(px_new)
        movegui(hFig, px_new);
    else
        movegui(hFig, 'center');
    end
catch
    movegui(hFig, 'center');
end
        
set(hFig, 'visible', 'on');
fliped = false;

MX = glist;
Z = zeros(length(glist), length(cL));
c = c(cidx);
for k = 1:length(cL)
    Z(:, k) = median(Yori(:, c == k), 2);
end
[Z] = gui.i_norm4heatmap(Z);

% figure;
% h2=heatmap(strrep(cL,'_','\_'),MX,Z);
% h2.Title = 'Marker Gene Heatmap';
% h2.XLabel = 'Group';
% h2.YLabel = 'Marker Gene';
% h2.Colormap = parula;
% h2.GridVisible = 'off';
% h2.CellLabelColor='none';
% h2.ColorLimits=[min(Z(:)), max(Z(:))];


    function i_flipxy(~, ~)
        %delete(h);
        fliped = ~fliped;
        if fliped
            h = imagesc(Y');
            set(gca, 'YTick', a-b);
            set(gca, 'YTickLabel', strrep(cL, '_', '\_'));
            %set(gca,'YTickLabelRotation',90);
            set(gca, 'XTick', 1:length(glist));
            set(gca, 'XTickLabel', glist);
            set(gca, 'XTickLabelRotation', 90);
            set(gca, 'TickLength', [0, 0]);
        else
            h = imagesc(Y);
            set(gca, 'XTick', a-b);
            set(gca, 'XTickLabel', strrep(cL, '_', '\_'));
            %set(gca,'XTickLabelRotation',0);
            set(gca, 'YTick', 1:length(glist));
            set(gca, 'YTickLabel', glist);
            set(gca, 'TickLength', [0, 0]);
        end
end

        function i_renamecat(~, ~)
            tg = gui.i_inputgenelist(string(cL), true);
            if isempty(tg), return; end
            if length(tg) == length(cL)
                set(gca, 'XTick', a-b);
                set(gca, 'XTickLabel', tg(:))
                cL = tg;
            else
                errordlg('Wrong input.');
            end
    end

            function i_resetcolor(~, ~)
                set(gca, 'FontSize', 10);
                colormap default
        end

                function i_summarymap(~, ~)
                    f = figure;
                    h = heatmap(strrep(cL, '_', '\_'), MX, Z);
                    h.Title = 'Marker Gene Heatmap';
                    h.XLabel = 'Group';
                    h.YLabel = 'Marker Gene';
                    h.Colormap = parula;
                    h.GridVisible = 'off';
                    h.CellLabelColor = 'none';
                    tb = uitoolbar('Parent', f);

                    t = array2table(Z, 'VariableNames', cL, 'RowNames', MX);

                    % writetable(t,'aaa.csv','WriteRowNames',true);
                    pkg.i_addbutton2fig(tb, 'off', {@i_exporttable, t}, 'export.gif', 'Save table...');

                    pkg.i_addbutton2fig(tb, 'on', {@gui.i_pickcolormap, c}, 'plotpicker-compass.gif', 'Pick new color map...');
                    pkg.i_addbutton2fig(tb, 'off', @gui.i_changefontsize, 'noun_font_size_591141.gif', 'ChangeFontSize');
                    pkg.i_addbutton2fig(tb, 'on', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
                    pkg.i_addbutton2fig(tb, 'on', @gui.i_invertcolor, 'plotpicker-comet.gif', 'Invert colors');
                    pkg.i_addbutton2fig(tb, 'off', @i_resetcolor, 'plotpicker-geobubble2.gif', 'Reset color map');
                    disp('https://software.broadinstitute.org/morpheus/')
            end

                    function i_exporttable(~, ~, T, needwait, defname)
                        if nargin < 5, defname = []; end
                        if nargin < 4, needwait = false; end
                        if ~isempty(defname)
                            [file, path] = uiputfile({'*.txt'; '*.*'}, 'Save as', defname);
                        else
                            [file, path] = uiputfile({'*.txt'; '*.*'}, 'Save as');
                        end
                        if isequal(file, 0) || isequal(path, 0)
                            return;
                        else
                            filename = fullfile(path, file);
                            try
                                writetable(T, filename, 'Delimiter', '\t', 'WriteRowNames', true);
                            catch
                                writematrix(T, filename, 'Delimiter', '\t');
                            end
                            pause(1);
                            if needwait
                                waitfor(helpdlg(sprintf('Result has been saved in %s', filename), ''));
                            else
                                helpdlg(sprintf('Result has been saved in %s', filename), '')
                            end
                        end
                end

                        function i_summarymapT(~, ~)
                            f = figure;
                            h = heatmap(MX, strrep(cL, '_', '\_'), Z.');
                            h.Title = 'Marker Gene Heatmap';
                            h.YLabel = 'Group';
                            h.XLabel = 'Marker Gene';
                            h.Colormap = parula;
                            h.GridVisible = 'off';
                            h.CellLabelColor = 'none';
                            %         s = struct(h);
                            %         s.XAxis.TickLabelRotation=45;
                            tb = uitoolbar('Parent', f);

                            t = array2table(Z.', 'VariableNames', MX, 'RowNames', cL);
                            % writetable(t,'aaa.csv','WriteRowNames',true);
                            pkg.i_addbutton2fig(tb, 'off', {@i_exporttable, t}, 'export.gif', 'Save table...');

                            pkg.i_addbutton2fig(tb, 'on', {@gui.i_pickcolormap, c}, 'plotpicker-compass.gif', 'Pick new color map...');
                            pkg.i_addbutton2fig(tb, 'off', @gui.i_changefontsize, 'noun_font_size_591141.gif', 'ChangeFontSize');
                            pkg.i_addbutton2fig(tb, 'on', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
                            pkg.i_addbutton2fig(tb, 'on', @gui.i_invertcolor, 'plotpicker-comet.gif', 'Invert colors');
                            pkg.i_addbutton2fig(tb, 'off', @i_resetcolor, 'plotpicker-geobubble2.gif', 'Reset color map');
                    end


end

