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

szgn = grpstats(c, c, @numel);
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

hx=gui.myFigure;
hFig=hx.FigHandle;
h = imagesc(Y);
set(gca, 'XTick', a-b);
set(gca, 'XTickLabel', strrep(cL, '_', '\_'));
%set(gca,'XTickLabelRotation',0);
set(gca, 'YTick', 1:length(glist));
set(gca, 'YTickLabel', glist);
set(gca, 'TickLength', [0, 0]);
% colormap(flipud(bone));
box on

szc = cumsum(szgn);
for kx = 1:length(szc)
    xline(szc(kx)+0.5, 'y-');
end


hx.addCustomButton('off', @i_renamecat, 'edit.jpg', 'Rename groups...');
hx.addCustomButton('off', @i_resetcolor, 'refresh_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Reset color map');
hx.addCustomButton('off', @i_flipxy, 'mat-wrap-text.jpg', 'Flip XY');
hx.addCustomButton('on', @i_summarymap, 'HDF_object01.gif', 'Summary map...');
hx.addCustomButton('off', @i_summarymapT, 'HDF_object02.gif', 'Summary map, transposed...');
hx.addCustomButton('on', @in_savetable, 'floppy-disk-arrow-in.jpg', 'Export data...');
hx.addCustomButton('off', @in_changenorm, 'mw-pickaxe-mining.jpg', 'Change normalization method...');
hx.addCustomButton('off', @i_dotplotx, 'icon-mat-blur-linear-10.gif', 'Dot plot...');

hx.show(parentfig);        
fliped = false;

MX = glist;

c = c(cidx);
Z = zeros(length(glist), length(cL));

% for k = 1:length(cL)
%     Z(:, k) = mean(Y(:, c == k), 2);
% end

%[Z] = gui.i_norm4heatmap(Z);

% figure;
% h2=heatmap(strrep(cL,'_','\_'),MX,Z);
% h2.Title = 'Marker Gene Heatmap';
% h2.XLabel = 'Group';
% h2.YLabel = 'Marker Gene';
% h2.Colormap = parula;
% h2.GridVisible = 'off';
% h2.CellLabelColor='none';
% h2.ColorLimits=[min(Z(:)), max(Z(:))];

    function in_changenorm(~, ~)
        [methodid, dim] = gui.i_selnormmethod;

        % [Xt] = gui.i_transformx(sce.X, true, 8, parentfig);
        % if isempty(Xt), return; end
        %  Yt = Xt(gidx, :);
        %  [~, cidx] = sort(c);
        %  Yori = Yt(:, cidx);

        [Y] = gui.i_norm4heatmap(Yori, dim, methodid);
        % Y = log1p(Y);
        delete(h);        
        h = imagesc(Y);
        set(gca, 'XTick', a-b);
        set(gca, 'XTickLabel', strrep(cL, '_', '\_'));
        set(gca, 'YTick', 1:length(glist));
        set(gca, 'YTickLabel', glist);
        set(gca, 'TickLength', [0, 0]);
        % colormap(flipud(bone));
        box on        
    end

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
            gui.myErrordlg(parentfig, 'Wrong input.');
        end
    end    

    function i_resetcolor(~, ~)
        set(gca, 'FontSize', 10);
        colormap default
    end        

    function i_summarymap(~, ~)
        for ky = 1:length(cL)
            Z(:, ky) = mean(Y(:, c == ky), 2);
        end
        
        hx1=gui.myFigure;
        h = heatmap(strrep(cL, '_', '\_'), MX, Z);
        h.Title = 'Marker Gene Heatmap';
        h.XLabel = 'Group';
        h.YLabel = 'Marker Gene';
        h.Colormap = parula;
        h.GridVisible = 'off';
        h.CellLabelColor = 'none';
        t = array2table(Z, 'VariableNames', cL, 'RowNames', MX);
        % writetable(t,'aaa.csv','WriteRowNames',true);
        hx1.addCustomButton('off', {@i_exporttable, t}, 'floppy-disk-arrow-in.jpg', 'Save table...');
        hx1.addCustomButton('off', @i_resetcolor, 'refresh_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Reset color map');
        % disp('https://software.broadinstitute.org/morpheus/')
        hx1.show(hFig);
    end

    function in_savetable(~, ~)
        labels = {'Save Y to variable named:', ...
            'Save glist to variable named:', ...
            'Save cL to variable named:'};
        vars = {'Y', 'g', 'cL'};
        values = {full(Y), glist, string(cL)};
        [~, ~] = export2wsdlg(labels, vars, values, ...
            'Save Data to Workspace');
    end

    function i_exporttable(~, ~, T, needwait, defname)
        if nargin < 5, defname = []; end
        if nargin < 4, needwait = false; end
        if ~isempty(defname)
            [file, path] = uiputfile({'*.txt'; '*.*'}, 'Save as', defname);
        else
            [file, path] = uiputfile({'*.txt'; '*.*'}, 'Save as');
        end
        if isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure'), figure(parentfig); end
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
                gui.myHelpdlg(parentfig, ...
                    sprintf('Result has been saved in %s', filename));
            end
        end
    end

    function i_summarymapT(~, ~)
            for ky = 1:length(cL)
                Z(:, ky) = mean(Y(:, c == ky), 2);
            end
        
        hx2=gui.myFigure;
        h = heatmap(MX, strrep(cL, '_', '\_'), Z.');
        h.Title = 'Marker Gene Heatmap';
        h.YLabel = 'Group';
        h.XLabel = 'Marker Gene';
        h.Colormap = parula;
        h.GridVisible = 'off';
        h.CellLabelColor = 'none';
        t = array2table(Z.', 'VariableNames', MX, 'RowNames', cL);
        %         s = struct(h);
        %         s.XAxis.TickLabelRotation=45;
        % writetable(t,'aaa.csv','WriteRowNames',true);
        hx2.addCustomButton('off', {@i_exporttable, t}, 'floppy-disk-arrow-in.jpg', 'Save table...');
        hx2.addCustomButton('off', {@gui.i_pickcolormap, c}, 'color-wheel.jpg', 'Pick new color map...');
        hx2.addCustomButton('off', @i_resetcolor, 'refresh_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Reset color map');
        hx2.show(hFig);
    end

    function i_dotplotx(~, ~)
        try
            gui.i_dotplot(sce.X, sce.g, c, cL, MX);
        catch ME
            gui.myErrordlg(parentfig, ME.message, ME.identifier);
        end
    end
    
end
