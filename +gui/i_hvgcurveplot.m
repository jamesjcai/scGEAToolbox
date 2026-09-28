function [xyz, xyz1, gsorted] = i_hvgcurveplot(X, g, dofit, showdata, ...
    parentfig, method)
% I_HVGCURVEPLOT 3-D mean/CV/dropout gene cloud with a fitted curve
%
%   METHOD picks the reference curve and, with it, the coordinates the
%   cloud is drawn in. The two cannot be chosen separately: a curve only
%   lies among the points it was fitted to.
%
%     "analytic" (default) - the closed-form gamma-Poisson curve of
%        SC_ANALYTICFIT, over library-size normalized coordinates. This
%        is the ranker SINGLECELLEXPERIMENT.EMBEDCELLS selects genes
%        with, so the picture explains the genes an embedding used.
%     "splinefit" - the smoothing spline of SC_SPLINEFIT over raw
%        counts, which is what this function drew before and what the
%        Splinefit option of GUI.CALLBACK_ENRICHRHVGS still asks for.
%
%   The two coordinate spaces are not interchangeable. On the bundled
%   8260-cell set the raw and normalized log mean differ by up to 0.978
%   over a 0.076-7.481 range, and the log CV by up to 0.484 over
%   0.472-3.271, so either curve drawn over the other's cloud sits off
%   the points it is meant to explain. Only dropout agrees, scaling not
%   being able to turn a zero into a nonzero.
%
% See also SC_ANALYTICFIT, SC_SPLINEFIT, SC_GENESTAT.

if nargin < 6 || isempty(method), method = "analytic"; end
if nargin < 5, parentfig = []; end
if nargin < 4, showdata = true; end
if nargin < 3, dofit = true; end
if nargin < 2 || isempty(g), g = string(1:size(X,1)); end
method = validatestring(method, ["analytic", "splinefit"], ...
    mfilename, 'method', 6);
% XYZ1 only: both DOFIT branches set XYZ, but a splinefit fit is the one
% thing that may never run, so its curve stays empty.
xyz1=[];

% Cloud and curve come from one fit, so they cannot fall out of step.
% SC_ANALYTICFIT returns gene names rather than a sorted matrix, so
% XSORTED is rebuilt by name; a repeated gene name takes its first row,
% which only decides which profile the data tip draws.
switch method
    case "analytic"
        [Tfit, xyz1, curve] = sc_analyticfit(X, g, SortIt=false);
        gsorted = Tfit.genes;
        lgu = Tfit.lgu; lgcv = Tfit.lgcv; dropr = Tfit.dropr;
        [~, loc] = ismember(gsorted, string(g(:)));
        Xsorted = X(loc, :);
    case "splinefit"
        [lgu, dropr, lgcv, gsorted, Xsorted] = sc_genestat(X, g);
end

% SC_GENESTAT hands back sparse vectors for a sparse X, which the fit
% and the scatter below both want dense.
lgu = full(lgu); lgcv = full(lgcv); dropr = full(dropr);

x = lgu;
y = lgcv;
z = dropr;

fw = gui.myWaitbar(parentfig);


hx=gui.myFigure(parentfig, true);
hFig=hx.FigHandle;
hFig.Position(3) = hFig.Position(3)*1.8;


hx.addCustomButton('off', @in_callback_HighlightTopHVGs, 'plotpicker-qqplot.gif', 'Highlight top HVGs');
hx.addCustomButton('off', {@in_HighlightSelectedGenes,2}, 'curve-array.jpg', 'Select HVG to show');
hx.addCustomButton('off', {@in_HighlightSelectedGenes,1}, 'checklist_rtl_18dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Highlight selected genes');
hx.addCustomButton('on', @in_callback_ExportTable, 'floppy-disk-arrow-in.jpg', 'Export HVG Table...');
hx.addCustomButton('off', @ExportGeneNames, 'bookmark-book.jpg', 'Export selected HVG gene names...');
hx.addCustomButton('off', @EnrichrHVGs, 'plotpicker-andrewsplot.gif', 'Enrichment analysis...');
hx.addCustomButton('off', @ChangeAlphaValue, 'Brightness-3--Streamline-Core.jpg', 'Change MarkerFaceAlpha value');

if showdata
    hAx1 = subplot(2,2,[1 3]);
    h = scatter3(hAx1, x, y, z, 'filled', 'MarkerFaceAlpha', .1);

    if ~isempty(gsorted)
        dt = datacursormode(hFig);
    else
        dt = [];
    end
end
    % grid on
    % box on
    % legend({'Genes','Spline fit'});
xlabel(hAx1,'Mean+1, log');
ylabel(hAx1,'CV+1, log');
zlabel(hAx1,'Dropout rate (% of zeros)');


if dofit
    switch method
        case "analytic"
            % SC_ANALYTICFIT has already scored every gene against its
            % own curve, deflating the ones below it exactly as the
            % splinefit branch does, so there is nothing to recompute.
            d = Tfit.d;

            % Draw the curve itself, not the per-gene foot points: in
            % gene order those zigzag, consecutive genes not being
            % consecutive along the curve. The grid spans the observed
            % means, so the line covers the cloud and no more.
            u = expm1(x);
            u = u(u > 0);
            mus = exp(linspace(log(min(u)), log(max(u)), 500)).';
            cxyz = [curve.x(mus), curve.y(mus), curve.z(mus)];
        case "splinefit"
            [~, ~, ~, xyz1] = sc_splinefit(Xsorted, gsorted);

            [~, d] = dsearchn(xyz1, [x, y, z]);
            fitmeanv=xyz1(:,1);
            d(x>max(fitmeanv))=d(x>max(fitmeanv))./100;
            d(x<min(fitmeanv))=d(x<min(fitmeanv))./10;
            d((y-xyz1(:, 2))<0)=d((y-xyz1(:, 2))<0)./100;

            % The spline is a curve through the gene cloud in gene
            % order, so its foot points already trace it.
            cxyz = xyz1;
    end

    hold(hAx1,'on');
    xyz = [x y z];
    plot3(hAx1, cxyz(:, 1), cxyz(:, 2), cxyz(:, 3), '-', 'linewidth', 4);

    [sortedd, hvgidx] = sort(d, 'descend');

    hvg=gsorted(hvgidx);
    lgu=lgu(hvgidx);
    lgcv=lgcv(hvgidx);
    dropr=dropr(hvgidx);

    gene=hvg(:);
    T=table(gene,sortedd,hvgidx,lgu,lgcv,dropr);

    disp('scGEAToolbox controls for the variance-mean relationship of gene')
    disp('expression. scGEAToolbox considers three sample statistics of each')
    disp('gene: expression mean, coefficient of variation, and dropout rate.')
    if method == "analytic"
        disp('It scores each gene by its distance to the closed-form')
        disp('gamma-Poisson curve those three statistics imply, fitted on')
        disp('library-size normalized counts. Genes with larger distances')
        disp('are ranked higher for feature selection.')
    else
        disp('After normalization, it fits a spline function based on')
        disp('piece-wise polynomials to model the relationship among the')
        disp('three statistics, and calculates the distance between each')
        disp('gene''s observed statistics and the fitted 3D spline curve.')
        disp('Genes with larger distances are ranked higher for feature')
        disp('selection.')
    end
else
    % No curve means no deviation to rank by, so gene order stands in for
    % the ranking. The subplot below indexes XSORTED with HVGIDX(1) and the
    % toolbar callbacks read T, both of which used to be unset here: a
    % DOFIT=false call threw "Unrecognized function or variable 'hvgidx'"
    % at the first SUBPLOT, which is why callers only ever passed true.
    xyz = [x y z];
    hvgidx = (1:numel(gsorted)).';
    hvg = gsorted;
    gene = hvg(:);
    sortedd = nan(numel(gsorted), 1);
    T = table(gene, sortedd, hvgidx, lgu, lgcv, dropr);
end


hAx2 = subplot(2,2,2);
x1=Xsorted(hvgidx(1),:);

sh = plot(hAx2, 1:length(x1), x1);
xlim(hAx2,[1 size(Xsorted,2)]);
title(hAx2, hvg(1));
[titxt] = gui.i_getsubtitle(x1);
subtitle(hAx2, titxt);
xlabel(hAx2,'Cell Index');
ylabel(hAx2,'Expression Level');

if showdata && ~isempty(dt)
    dt.UpdateFcn = {@in_myupdatefcn3, gsorted};
end
gui.myWaitbar(parentfig, fw);
hx.show(parentfig);


function ChangeAlphaValue(~, ~)
        if h.MarkerFaceAlpha <= 0.05
            h.MarkerFaceAlpha = 1;
        else
            h.MarkerFaceAlpha = h.MarkerFaceAlpha - 0.1;
        end
    end

function in_callback_HighlightTopHVGs(~, ~)
        idx = zeros(1, length(hvgidx));
        h.BrushData = idx;

        k = gui.i_inputnumk(200, 1, 2000, [], hFig);
        if isempty(k), return; end
        idx(hvgidx(1:k)) = 1;
        h.BrushData = idx;
    end

function in_callback_ExportTable(~, ~)
        gui.i_exporttable(T, true, 'Thvg', 'HVGTable', [], [], hFig);
    end

function in_HighlightSelectedGenes(~,~,typeid)
        if nargin<3, typeid = 1; end

           switch typeid
               case 1
                    % Myc, Oct3/4, Sox2, Klf4
                    [glist] = gui.i_selectngenes(SingleCellExperiment(Xsorted,gsorted),...
                        intersect(upper(gsorted),["MYC", "POU5F1", "SOX2", "KLF4"]), parentfig);
               case 2
                    % GLISTALL, not GSORTED: these are nested functions
                    % sharing one workspace, so assigning to GSORTED
                    % here replaced the scatter's gene order with the
                    % table's rank order, and every later ISMEMBER
                    % against it brushed the wrong points.
                    glistall = T.(T.Properties.VariableNames{1});

                   if gui.i_isuifig(parentfig)
                        [indx2, tf2] = gui.myListdlg(hFig, glistall, 'Select genes:', [], true);
                    else
                        [indx2, tf2] = listdlg('PromptString', ...
                            'Select genes:', ...
                            'SelectionMode', 'multiple', ...
                            'ListString', glistall, ...
                            'ListSize', [220, 300]);
                    end


                    if tf2 == 1
                        glist = glistall(indx2);
                    else
                        return;
                    end
           end

        if ~isempty(glist)
            [yes,idx]=ismember(glist,gsorted);
            idx=idx(yes);

            for k=1:length(idx)
                dt = datatip(h,'DataIndex',idx(k));
            end
        end
    end

function ExportGeneNames(~, ~)
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            gui.myWarndlg(hFig, "No gene is selected.");
            return;
        end
        fprintf('%s selected.\n', pkg.i_plural(sum(ptsSelected), 'gene'));

        gselected=gsorted(ptsSelected);
        [yes,idx]=ismember(gselected,T.gene);
        Tx=T(idx,:);
        Tx=sortrows(Tx,1,'descend');
        if ~all(yes), error('Running time error.'); end
        tgenes=Tx.gene;

        labels = {'Save gene names to variable:'};
        vars = {'g'};
        values = {tgenes};
        export2wsdlg(labels, vars, values, ...
            'Save Data to Workspace');
    end

function EnrichrHVGs(~, ~)
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            gui.myWarndlg(hFig, "No gene is selected.");
            return;
        end
        fprintf('%s selected.\n', pkg.i_plural(sum(ptsSelected), 'gene'));

        gselected=gsorted(ptsSelected);
        [yes,idx]=ismember(gselected,T.gene);
        Tx=T(idx,:);
        Tx=sortrows(Tx,1,'descend');
        if ~all(yes), error('Running time error.'); end
        tgenes=Tx.gene;
        gui.i_enrichtest(tgenes, gsorted, numel(tgenes));
    end

function txt = in_myupdatefcn3(src, event_obj, g)
        if isequal(get(src, 'Parent'), hAx1)
            idx = event_obj.DataIndex;
            txt = {g(idx)};
            x1 = Xsorted(idx, :);

            if ~isempty(sh) && pkg.i_isvalid(sh)
                delete(sh);
            end
            sh = plot(hAx2, 1:length(x1), x1, 'marker', 'none');
            xlim(hAx2,[1 size(Xsorted,2)]);
            title(hAx2, g(idx));
            [titxt] = gui.i_getsubtitle(x1);
            subtitle(hAx2, titxt);
            xlabel(hAx2,'Cell Index');
            ylabel(hAx2,'Expression Level');
        else
            txt = num2str(event_obj.Position(2));
        end
    end

end
