function callback_FindAllMarkers(src, ~)

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

answer = questdlg('Select Method', ...
    '', 'Find All Markers', 'Marker Gene Heatmap', ...
    'Find All Markers');
switch answer
    case 'Find All Markers'

    case 'Marker Gene Heatmap'
        i_MarkerGeneHeatmap(src);
        return;
    otherwise
        return;
end

[thisc, ~] = gui.i_select1class(sce);
if isempty(thisc), return; end
if numel(unique(thisc))==1
    warndlg("All cells are in the same group.",'');
    return;
end

[T] = pkg.e_findallmarkers(sce.X, sce.g, thisc, [], [], [], true);
if ~isempty(T)
    needwait = true;
    [answer, filename] = gui.i_exporttable(T, needwait, 'Tallmarkers', ...
        'AllMarkersTable');
            % "Tcellattrib","CellAttribTable"
            % "Tviolindata","ViolinPlotTable"
            % "Tcrosstabul","CrosstabulTable"
            % "Tcellsignmt","CellSignatTable"
            % "Tdpgenesres","DPGenesResTable"
            % "Tallmarkers","AllMarkersTable"
    if ~isempty(answer)
        disp(filename);
        helpdlg(sprintf('All Markers Table saved.'), '');
    end
else
    helpdlg('No results.', '');
end
end


function i_MarkerGeneHeatmap(src, ~, sce)

if nargin < 3
    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
end
% unique(sce.c_cluster_id)

[thisc, ~] = gui.i_select1class(sce);
if isempty(thisc), return; end
if numel(unique(thisc))==1
    warndlg("All cells are in the same group.",'');
    return;
end
[c, cL, noanswer] = gui.i_reordergroups(thisc);
if noanswer, return; end

answer = questdlg('Generate marker gene heatmap', ...
    'Select Method', 'Method 1 (DE 🐇)', 'Method 2 (scGeneFit 🐢)', ...
    'Method 3 (LASSO 🐢🐢)', 'Method 1 (DE 🐇)');
switch answer
    case 'Method 1 (DE 🐇)'
        methodid = 1;
    case 'Method 2 (scGeneFit 🐢)'
        methodid = 2;
    case 'Method 3 (LASSO 🐢🐢)'
        methodid = 3;
    otherwise
        return;
end

fw = gui.gui_waitbar;
try
    [markerlist] = sc_pickmarkers(sce.X, sce.g, c, 10, methodid);
catch ME
    gui.gui_waitbar(fw, true);
    errordlg(ME.message);
    return;
end

M = cell(numel(cL), 2);
for k = 1:numel(cL)
    %cLk=matlab.lang.makeValidName(cL{k});
    cLk = cL{k};
    M{k, 1} = cLk;
    M{k, 2} = markerlist{k};
end

X = [];
szcl = [];
idcl = [];
for k = 1:length(cL)
    i = c == k;
    X = [X, sce.X(:, i)];
    szcl = [szcl, sum(i)];
    idcl = [idcl; c(i)];
end
X = sc_norm(X);
X = log(X+1);

% ===========
Y = [];
idgn = [];
szgn = [];
Z = [];
% subi=1:10:size(X,2);
MX = [];
for k = 1:numel(cL)
    markerlist = M{k, 2}(1:end);
    MX = [MX; markerlist];
    [~, idx_g] = ismember(upper(markerlist), upper(sce.g));
    Y = [Y; X(idx_g, :)];
    idgn = [idgn; k * ones(length(markerlist), 1)];
    szgn = [szgn, length(markerlist)];
end
Yraw = Y;
[Y] = gui.i_norm4heatmap(Y);

Z = [];
Zfd = [];
for kx = 1:numel(cL)
    y = Y(idgn == kx, :);
    yraw = Yraw(idgn == kx, :);
    z = [];
    zfd = [];
    for kk = 1:numel(cL)
        z = [z, mean(y(:, idcl == kk), 2)];
        zfd = [zfd, log2(mean(yraw(:, idcl == kk), 2)./mean(yraw(:, idcl ~= kk), 2))];
    end
    %z1=grpstats(y.',idcl,@mean)';
    %assert(isequal(z,z1));
    Z = [Z; z];
    Zfd = [Zfd; zfd];
end
gui.gui_waitbar(fw);


% ======= customized heatmap - start
hFig = figure('Visible', 'off');
himg = imagesc(Y);
szc = cumsum(szgn);
for kx = 1:max(idcl) - 1
    xline(sum(idcl < kx+1)+0.5, 'r-');
    yline(szc(kx)+0.5, 'r-');
end

set(gca, 'YTick', 1:size(Y, 1));
a = zeros(1, max(idcl));
b = zeros(1, max(idcl));
for kx = 1:max(idcl)
    a(kx) = sum(idcl <= kx);
    b(kx) = round(sum(idcl == kx)./2);
end

set(gca, 'XTick', a-b);
% set(gca,'XTickLabel',strrep(M(:,1),'_','\_'));
set(gca, 'XTickLabel', M(:, 1));
%set(gca,'XTickLabelRotation',0);
set(gca, 'YTick', 1:length(MX));
set(gca, 'YTickLabel', MX);
set(gca, 'TickLength', [0, 0])
% ======= customized heatmap - end

% tb1 = uitoolbar(hFig);
tb1 = findall(hFig, 'Tag', 'FigureToolBar'); % get the figure's toolbar handle
uipushtool(tb1, 'Separator', 'off');

pkg.i_addbutton2fig(tb1, 'off', {@i_saveM, M}, 'greencircleicon.gif', 'Save marker gene map...');
pkg.i_addbutton2fig(tb1, 'off', @i_flipxy, 'xplotpicker-geobubble2.gif', 'Flip XY');
pkg.i_addbutton2fig(tb1, 'off', {@i_summarymap, Z}, 'HDF_object01.gif', 'Summary map...');
pkg.i_addbutton2fig(tb1, 'off', {@i_summarymap, Z.'}, 'HDF_object02.gif', 'Summary map, transposed...');
pkg.i_addbutton2fig(tb1, 'off', @i_dotplotx, 'HDF_object03.gif', 'Dot plot...');
pkg.i_addbutton2fig(tb1, 'on', {@gui.i_pickcolormap, c}, 'plotpicker-compass.gif', 'Pick new color map...');
pkg.i_addbutton2fig(tb1, 'off', @gui.i_changefontsize, 'noun_font_size_591141.gif', 'ChangeFontSize');
pkg.i_addbutton2fig(tb1, 'on', @i_renamecat, 'guideicon.gif', 'Rename groups...');
pkg.i_addbutton2fig(tb1, 'on', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
pkg.i_addbutton2fig(tb1, 'on', {@i_savetable, Y}, 'export.gif', 'Export data...');

pkg.i_addbutton2fig(tb1, 'on', @gui.i_invertcolor, 'plotpicker-comet.gif', 'Invert colors');
pkg.i_addbutton2fig(tb1, 'off', @i_resetcolor, 'plotpicker-geobubble2.gif', 'Reset color map');
pkg.i_addbutton2fig(tb1, 'on', {@gui.i_resizewin, hFig}, 'HDF_pointx.gif', 'Resize Plot Window');

gui.i_movegui2parent(hFig, FigureHandle);

set(hFig, 'visible', 'on');
fliped = false;

    function i_flipxy(~, ~)
        fliped = ~fliped;
        if fliped
            himg = imagesc(Y');
            set(gca, 'YTick', a-b);
            set(gca, 'YTickLabel', M(:, 1));
            %set(gca,'YTickLabelRotation',90);
            set(gca, 'XTick', 1:length(MX));
            set(gca, 'XTickLabel', MX);
            set(gca, 'XTickLabelRotation', 90);
            set(gca, 'TickLength', [0, 0]);
        else
            himg = imagesc(Y);
            %set(gca,'YTick',1:size(Y,1));
            set(gca, 'XTick', a-b);
            set(gca, 'XTickLabel', M(:, 1));
            set(gca, 'YTick', 1:length(MX));
            set(gca, 'YTickLabel', MX);
            set(gca, 'TickLength', [0, 0])

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

                
    function i_invertcolor(~, ~)
        cm = colormap();
        colormap(flipud(cm));
    end


                    
    function i_dotplotx(~, ~)
        try
            gui.i_dotplot(sce.X, sce.g, c, cL, MX);
        catch ME
            % if exist('f', 'var') && ishandle(f)
            %     close(f);
            % end
            errordlg(ME.message);
        end
    end


    function i_summarymap(~, ~, thisZ)
        f = figure;
        if length(cL) == size(thisZ, 1)
            h = heatmap(MX, cL, thisZ);
        else
            h = heatmap(cL, MX, thisZ);
        end
        h.Title = 'Marker Gene Heatmap';
        h.XLabel = 'Group';
        h.YLabel = 'Marker Gene';
        h.Colormap = parula;
        h.GridVisible = 'off';
        h.CellLabelColor = 'none';
        tb = uitoolbar('Parent', f);
        pkg.i_addbutton2fig(tb, 'on', {@gui.i_pickcolormap, c}, 'plotpicker-compass.gif', 'Pick new color map...');
        pkg.i_addbutton2fig(tb, 'off', @gui.i_changefontsize, 'noun_font_size_591141.gif', 'ChangeFontSize');
        pkg.i_addbutton2fig(tb, 'on', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
        pkg.i_addbutton2fig(tb, 'on', @i_invertcolor, 'plotpicker-comet.gif', 'Invert colors');
        pkg.i_addbutton2fig(tb, 'off', @i_resetcolor, 'plotpicker-geobubble2.gif', 'Reset color map');
        pkg.i_addbutton2fig(tb, 'on', {@gui.i_resizewin, f}, 'HDF_pointx.gif', 'Resize Plot Window');
        if length(cL) ~= size(thisZ, 1)
            pkg.i_addbutton2fig(tb, 'on', {@i_savetable, thisZ}, 'export.gif', 'Export data...');
        end
    end


    function i_saveM(~, ~, M)
        if ~(ismcc || isdeployed)
            labels = {'Save marker gene map M to variable named:', ...
                'Save marker gene list G to variable named:'};
            vars = {'M', 'G'};
            G = cat(1, M{:, 2});
            values = {M, G};
            export2wsdlg(labels, vars, values);
        else
            errordlg('This function is not available for standalone application. Run scgeatool.m in MATLAB to use this function.');
        end        
    end


    function i_savetable(~, ~, c)
        answer = questdlg('Export & save data to:', '', ...
            'Workspace', 'TXT/CSV file', 'Excel file', 'Workspace');
        if ~isempty(answer)
            T = table();
            T.genes = cat(1, M{:, 2});
            assert(isequal(T.genes, MX))

            T = [T, array2table(c)];
            switch answer
                case 'Workspace'
                    labels = {'Save T to variable named:'};
                    vars = {'T'};
                    values = {T};
                    [~, ~] = export2wsdlg(labels, vars, values, ...
                        'Save Data to Workspace');
                case 'TXT/CSV file'
                    [file, path] = uiputfile({'*.csv'; '*.*'}, 'Save as');
                    if isequal(file, 0) || isequal(path, 0)
                        return;
                    else
                        fw = gui.gui_waitbar;
                        filename = fullfile(path, file);
                        writetable(T, filename, 'FileType', 'text');
                        OKPressed = true;
                        gui.gui_waitbar(fw);
                    end
                case 'Excel file'

                    [file, path] = uiputfile({'*.xlsx'; '*.*'}, 'Save as');
                    if isequal(file, 0) || isequal(path, 0)
                        return;
                    else
                        fw = gui.gui_waitbar;
                        filename = fullfile(path, file);
                        writetable(T, filename, 'FileType', 'spreadsheet');
                        OKPressed = true;
                        gui.gui_waitbar(fw);
                    end
            end
        end
    end

end