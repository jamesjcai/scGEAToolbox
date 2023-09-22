function sc_scatter3genes(X, g, dofit, showdata)
%Scatter3 plot for genes

if nargin < 4, showdata = true; end
if nargin < 3, dofit = true; end
if nargin < 2, g = []; end
[lgu, dropr, lgcv, g] = sc_genestat(X, g);
x = lgu;
y = lgcv;
z = dropr;
if showdata
    FigureHandle = figure;
    hAx = axes('Parent', FigureHandle);
    UitoolbarHandle = uitoolbar('Parent', FigureHandle);
    set(UitoolbarHandle, 'Tag', 'FigureToolBar', ...
        'HandleVisibility', 'off', 'Visible', 'on');


    pkg.i_addbutton2fig(UitoolbarHandle, 'off', @HighlightGenes, 'plotpicker-qqplot.gif', 'Highlight top HVGs');
    pkg.i_addbutton2fig(UitoolbarHandle, 'off', @ExportGeneNames, 'export.gif', 'Export selected HVG gene names...');
    pkg.i_addbutton2fig(UitoolbarHandle, 'off', @EnrichrHVGs, 'plotpicker-andrewsplot.gif', 'Enrichment analysis...');
    pkg.i_addbutton2fig(UitoolbarHandle, 'off', @ChangeAlphaValue, 'xplotpicker-andrewsplot.gif', 'Change MarkerFaceAlpha value');

    gui.add_3dcamera(UitoolbarHandle, 'HVGs');
    pkg.i_addbutton2fig(UitoolbarHandle, 'off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');

    %h=scatter3(hAx,x,y,z);  % 'filled','MarkerFaceAlpha',.5);
    h = scatter3(hAx, x, y, z, 'filled', 'MarkerFaceAlpha', .1);

    if ~isempty(g)
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn1, g};
    end
end

%grid on
%box on
%legend({'Genes','Spline fit'});
xlabel('Mean, log');
ylabel('CV, log');
zlabel('Dropout rate (% of zeros)');


% [xData, yData, zData] = prepareSurfaceData(x,y,z);
% xyz=[xData yData zData]';
if dofit

    try
        [~, xyz1] = sc_splinefit(X, g);
    catch ME
        rethrow(ME);
    end

    %     xyz=[x y z]';
    %     % xyz=sortrows([x y z],[1 2])';
    %     pieces = 15;
    %     s = cumsum([0;sqrt(diff(x(:)).^2 + diff(y(:)).^2 + diff(z(:)).^2)]);
    %     pp1 = splinefit(s,xyz,pieces,0.75);
    %     xyz1 = ppval(pp1,s);
    hold on
    plot3(xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), '-', 'linewidth', 4);
    % scatter3(xyz1(:,1),xyz1(:,2),xyz1(:,3)); %,'MarkerEdgeAlpha',.8);

    [~, d] = dsearchn(xyz1, [x, y, z]);
    [~, hvgidx] = sort(d, 'descend');

    %g(idx20)

    disp('scGEAToolbox controls for the variance-mean relationship of gene')
        disp('expression. scGEAToolbox considers three sample statistics of each')
        disp('gene: expression mean⁠, coefficient of variation⁠, and dropout rate⁠.')
        disp('After normalization, it fits a spline function based on piece-wise')
            disp('polynomials to model the relationship among the three statistics, ')
            disp('and calculates the distance between each geneis observed statistics')
            disp('to the fitted 3D spline surface. Genes with larger distances are ')
            disp('ranked higher for feature selection.')

            end

                function ChangeAlphaValue(~, ~)
                    if h.MarkerFaceAlpha <= 0.05
                        h.MarkerFaceAlpha = 1;
                    else
                        h.MarkerFaceAlpha = h.MarkerFaceAlpha - 0.1;
                    end
            end

                    function HighlightGenes(~, ~)
                        %h.MarkerIndices=idx20;
                        k = gui.i_inputnumk(200, 10, 2000);
                        if isempty(k), return; end
                        idx = zeros(1, length(hvgidx));
                        idx(hvgidx(1:k)) = 1;
                        h.BrushData = idx;
                        % datatip(h, 'DataIndex', idx20);
                        %h2=scatter3(x(idx20),y(idx20),z(idx20),'rx');  % 'filled','MarkerFaceAlpha',.5);
                end

                        function ExportGeneNames(~, ~)
                            ptsSelected = logical(h.BrushData.');
                            if ~any(ptsSelected)
                                warndlg("No gene is selected.");
                                return;
                            end
                            fprintf('%d genes are selected.\n', sum(ptsSelected));

                            labels = {'Save gene names to variable:'};
                            vars = {'g'};
                            values = {g(ptsSelected)};
                            export2wsdlg(labels, vars, values, ...
                                'Save Data to Workspace');
                    end

                            function EnrichrHVGs(~, ~)
                                ptsSelected = logical(h.BrushData.');
                                if ~any(ptsSelected)
                                    warndlg("No gene is selected.");
                                    return;
                                end
                                fprintf('%d genes are selected.\n', sum(ptsSelected));

                                tgenes = g(ptsSelected);
                                answer = gui.timeoutdlg(@(x){questdlg('Which analysis?', '', ...
                                    'Enrichr', 'GOrilla', 'Enrichr+GOrilla', 'Enrichr')}, 15);
                                if isempty(answer), return; end
                                switch answer
                                    case 'Enrichr'
                                        run.web_Enrichr(tgenes, length(tgenes));
                                    case 'GOrilla'
                                        run.web_GOrilla(tgenes);
                                    case 'Enrichr+GOrilla'
                                        run.web_Enrichr(tgenes, length(tgenes));
                                        run.web_GOrilla(tgenes);
                                    otherwise
                                        return;
                                end

                        end

                        end


                        function txt = i_myupdatefcn1(~, event_obj, g)
                            % Customizes text of data tips
                            % pos = event_obj.Position;
                            idx = event_obj.DataIndex;
                            % i_plotsiglegene(idx,g);
                            txt = {g(idx)};
                        end