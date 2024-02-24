function sc_scatter3genes(X, g, dofit, showdata)
%Scatter3 plot for genes

if nargin < 4, showdata = true; end
if nargin < 3, dofit = true; end
if nargin < 2 || isempty(g), g = string(1:size(X,1)); end

[lgu, dropr, lgcv, g, X] = sc_genestat(X, g);

x = lgu;
y = lgcv;
z = dropr;
if showdata
    FigureHandle = figure;
    hAx = axes('Parent', FigureHandle);
    UitoolbarHandle = uitoolbar('Parent', FigureHandle);
    set(UitoolbarHandle, 'Tag', 'FigureToolBar', ...
        'HandleVisibility', 'off', 'Visible', 'on');

    pkg.i_addbutton2fig(UitoolbarHandle, 'off', @in_ShowProfile, 'plotpicker-qqplotx.gif', 'Show Profile of Genes');

    pkg.i_addbutton2fig(UitoolbarHandle, 'off', @in_HighlightGenes, 'plotpicker-qqplot.gif', 'Highlight top HVGs');
    pkg.i_addbutton2fig(UitoolbarHandle, 'off', @in_HighlightSelectedGenes, 'xplotpicker-qqplot.gif', 'Highlight selected genes');
    pkg.i_addbutton2fig(UitoolbarHandle, 'off', @ExportGeneNames, 'export.gif', 'Export selected HVG gene names...');
    pkg.i_addbutton2fig(UitoolbarHandle, 'off', @ExportTable, 'xexport.gif', 'Export HVG Table...');
    pkg.i_addbutton2fig(UitoolbarHandle, 'off', @EnrichrHVGs, 'plotpicker-andrewsplot.gif', 'Enrichment analysis...');
    pkg.i_addbutton2fig(UitoolbarHandle, 'off', @ChangeAlphaValue, 'xplotpicker-andrewsplot.gif', 'Change MarkerFaceAlpha value');

    gui.add_3dcamera(UitoolbarHandle, 'HVGs');
    pkg.i_addbutton2fig(UitoolbarHandle, 'off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');

    %h=scatter3(hAx,x,y,z);  % 'filled','MarkerFaceAlpha',.5);
    h = scatter3(hAx, x, y, z, 'filled', 'MarkerFaceAlpha', .1);

    if ~isempty(g)
        dt = datacursormode(FigureHandle);
        dt.UpdateFcn = {@i_myupdatefcn1, g};
    else
        dt = [];
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
        [~, ~, ~, xyz1] = sc_splinefit(X, g);
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


    fitmeanv=xyz1(:,1);
    d(x>max(fitmeanv))=d(x>max(fitmeanv))./100;
    d(x<min(fitmeanv))=d(x<min(fitmeanv))./10;
    d((y-xyz1(:, 2))<0)=d((y-xyz1(:, 2))<0)./100;

    [sortedd, hvgidx] = sort(d, 'descend');

    hvg=g(hvgidx);
    lgu=lgu(hvgidx);
    lgcv=lgcv(hvgidx);
    dropr=dropr(hvgidx);    

    T=table(sortedd,hvgidx,hvg,lgu,lgcv,dropr);
    %assignin("base","T",T);
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

    function in_ShowProfile(~, ~)
        if ~isempty(dt)
           dt.UpdateFcn = {@i_myupdatefcn3, g, X};
        end
    end

    function ChangeAlphaValue(~, ~)
        if h.MarkerFaceAlpha <= 0.05
            h.MarkerFaceAlpha = 1;
        else
            h.MarkerFaceAlpha = h.MarkerFaceAlpha - 0.1;
        end
    end            

    function in_HighlightGenes(~, ~)
        %h.MarkerIndices=idx20;
        idx = zeros(1, length(hvgidx));
        h.BrushData = idx;

        k = gui.i_inputnumk(200, 1, 2000);
        if isempty(k), return; end
        

        idx(hvgidx(1:k)) = 1;
        h.BrushData = idx;
        % datatip(h, 'DataIndex', idx20);
        %h2=scatter3(x(idx20),y(idx20),z(idx20),'rx');  % 'filled','MarkerFaceAlpha',.5);
    end

    function ExportTable(~, ~)                
        gui.i_exporttable(T, true, 'Tsplinefitg', 'SplinefitGTable');
    end

    function in_HighlightSelectedGenes(~,~)        
        %Myc, Oct3/4, Sox2, Klf4
        [glist] = gui.i_selectngenes(SingleCellExperiment(X,g),...
            intersect(upper(g),["MYC", "POU5F1", "SOX2", "KLF4"]));
        if ~isempty(glist)            
            [yes,idx]=ismember(glist,g);
            idx=idx(yes);

            %idx=[idx(:); find(nearestidx==1)];
            % idv = zeros(1, length(hvgidx));
            % idv(idx)=1;
            % h.BrushData = idv;
            for k=1:length(idx)
                dt = datatip(h,'DataIndex',idx(k));
            end
        end
    end

    function ExportGeneNames(~, ~)
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warndlg("No gene is selected.");
            return;
        end
        fprintf('%d genes are selected.\n', sum(ptsSelected));

        gselected=g(ptsSelected);
        [yes,idx]=ismember(gselected,T.hvg);
        Tx=T(idx,:);
        Tx=sortrows(Tx,1,'descend');        
        if ~all(yes), error('Running time error.'); end
        tgenes=Tx.hvg;

        labels = {'Save gene names to variable:'};
        vars = {'g'};
        values = {tgenes};
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

        gselected=g(ptsSelected);
        [yes,idx]=ismember(gselected,T.hvg);
        Tx=T(idx,:);
        Tx=sortrows(Tx,1,'descend');        
        if ~all(yes), error('Running time error.'); end
        tgenes=Tx.hvg;

        gui.i_enrichtest(tgenes, g, numel(tgenes));

        % answer = gui.timeoutdlg(@(x){questdlg('Which analysis?', '', ...
        %     'Enrichr', 'GOrilla', 'Enrichr+GOrilla', 'Enrichr')}, 15);
        % if isempty(answer), return; end
        % switch answer
        %     case 'Enrichr'
        %         run.web_Enrichr(tgenes, length(tgenes));
        %     case 'GOrilla'
        %         run.web_GOrilla(tgenes);
        %     case 'Enrichr+GOrilla'
        %         run.web_Enrichr(tgenes, length(tgenes));
        %         run.web_GOrilla(tgenes);
        %     otherwise
        %         return;
        % end
    end        
end
