function [T, Xsorted, gsorted] = sc_hvg(X, g, sortit, plotit, ...
    normit, ignorehigh, ignorelow)
% Identify HVGs
% HVGs selection - This method uses the CV^2 on normalized count data to
% analyze biological variation.
%
% REF: https://www.nature.com/articles/nmeth.2645
% Input X: Brennecke et al. (2013) uses DESeq's method for normalization.
%
% USAGE:
% >> [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% >> [T]=sc_hvg(X,genelist);
%
% See also: SC_SPLINEFIT, SC_VEG

if nargin < 2 || isempty(g)
    g = strcat("G", string(1:size(X, 1)))';
end
if nargin < 3, sortit = true; end
if nargin < 4, plotit = false; end
if nargin < 5, normit = true; end
if nargin < 6, ignorehigh = true; end
if nargin < 7, ignorelow = true; end

%if nargout > 1, 
    Xori = X; 
%end


dropr = 1 - sum(X > 0, 2) ./ size(X, 2);
if normit
    %[X]=pkg.norm_deseq(X);
    [X] = pkg.norm_libsize(X);
end
if any(isnan(X(:)))
    u = mean(X, 2, 'omitnan');
    vx = var(X, 0, 2, 'omitnan');
    cv2 = vx ./ u.^2;
else
    % issparse(X)
    u = mean(X, 2);
    vx = var(X, 0, 2);
    % vx=sum(abs(X-mean(X,2)).^2,2)./(size(X,2)-1);
    cv2 = vx ./ u.^2;    
end

if issparse(u), u = full(u); end
if issparse(vx), vx = full(vx); end
if issparse(cv2), cv2 = full(cv2); end

xi = 1 ./ u;
yi = cv2;

removedidx1 = false(length(xi),1);
removedidx2 = false(length(xi),1);
if ignorelow
    % removing genes with the dropoff rate >= 0.95
    idx1 = dropr >= 0.95;
    removedidx1(idx1) = true;
end

if ignorehigh
    removedidx2(u > quantile(u, 0.99)) = true;
end

xi(removedidx1|removedidx2)=[];
yi(removedidx1|removedidx2)=[];

m = size(X, 2);
df = m - 1;

% b=glmfit(xi,yi,'gamma','link','identity');
% cv2fit=glmval(b,1./u,'identity');    % OR cv2fit=b(2)./u+b(1);
% if issparse(xi), xi=full(xi); end
% if issparse(yi), yi=full(yi); end
mdl = fitglm(xi, yi, 'linear', 'Distribution', 'gamma', 'link', 'identity');
cv2fit = mdl.predict(1./u);
b = mdl.Coefficients.Estimate;


methodid = 1;
switch methodid
    case 1
        % this code follows: https://github.com/tallulandrews/M3Drop/blob/master/R/Brennecke_implementation.R
        % and pages 7-9 https://media.nature.com/original/nature-assets/nmeth/journal/v10/n11/extref/nmeth.2645-S2.pdf
        minBiolDisp = 0.5.^2;
        cv2th = b(1) + minBiolDisp + b(1) * minBiolDisp;
        testDenom = (u * b(2) + u.^2 * cv2th) / (1 + cv2th / m);
        fitratio = vx ./ testDenom;
    case 2
        % this code follows https://github.com/MarioniLab/MNN2017/blob/a202f960f165816f22dec3b62ce1c7549b3ba8c1/Pancreas/findHighlyVariableGenes.R
        % cv2fit=b(2)./u+b(1);
        fitratio = cv2 ./ cv2fit;
end

pval = chi2cdf(fitratio*df, df, 'upper');
% OR 1-chi2cdf(fitratio*df,df);
residualcv2 = log(fitratio); % log(cv2)-log(cv2fit);

% fdr=mafdr(pval,'BHFDR',true);
[~, ~, ~, fdr] = pkg.fdr_bh(pval);

T = table(g, u, cv2, residualcv2, dropr, fitratio, pval, fdr, removedidx1, removedidx2);


T.Properties.VariableNames(1) = {'genes'};

%i = ~isnan(cv2);
%T = T(i, :);

if sortit
    T.fitratio(T.dropr > (1 - 0.05)) = 0; % ignore genes with dropout rate > 0.95
    % disp('NOTE: Genes with dropout rate > 0.95 are excluded.');
    [T, hvgidx] = sortrows(T, 'fitratio', 'descend');
    %if nargout > 1
        Xsorted = Xori(hvgidx, :);
    %end
    %if nargout > 2
        gsorted = T.genes;
    %end
else
    %if nargout > 1
        error('SORTIT was not required.');
    %end
end
% T=T(removedidx,:);
if plotit
    [~,top100idx]=maxk(fitratio,100);
    %    figure;

    hFig = figure;
    hFig.Position(3)=hFig.Position(3)*1.8;
    %hAx = axes('Parent', FigureHandle);
    hAx1 = subplot(2,2,[1 3]);
    % tb = uitoolbar('Parent', FigureHandle);
    % set(tb, 'Tag', 'FigureToolBar', ...
    %     'HandleVisibility', 'off', 'Visible', 'on');
    tb = findall(hFig, 'Tag', 'FigureToolBar'); % get the figure's toolbar handle
uipushtool(tb, 'Separator', 'off');
    pkg.i_addbutton2fig(tb, 'off', @HighlightGenes, 'plotpicker-qqplot.gif', 'Highlight top HVGs');
    pkg.i_addbutton2fig(tb, 'off', @in_HighlightSelectedGenes, 'xplotpicker-qqplot.gif', 'Highlight selected genes');
    pkg.i_addbutton2fig(tb, 'off', @ExportGeneNames, 'export.gif', 'Export Selected HVG gene names...');
    pkg.i_addbutton2fig(tb, 'off', @ExportTable, 'xexport.gif', 'Export HVG Table...');
    pkg.i_addbutton2fig(tb, 'off', @EnrichrHVGs, 'plotpicker-andrewsplot.gif', 'Enrichment analysis...');
    pkg.i_addbutton2fig(tb, 'off', @ChangeAlphaValue, 'xplotpicker-andrewsplot.gif', 'Change MarkerFaceAlpha value');
    pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
    h = scatter(hAx1, log(u), log(cv2), 'filled', 'MarkerFaceAlpha', .1);
    hold on
    % scatter(log(u(top100idx)),log(cv2(top100idx)),'x');
    plot(hAx1, log(u), log(cv2fit), '.', 'markersize', 10);
    %plot(hAx1, log(u(removedidx1)), log(cv2(removedidx1)), 'xr', 'markersize', 10);
    %plot(hAx1, log(u(removedidx2)), log(cv2(removedidx2)), '+r', 'markersize', 10);
    %plot(hAx1, log(u(top100idx)), log(cv2(top100idx)), '^k', 'markersize', 10);

    %[~,i]=sort(fitratio,'descend');
    %xi=u(i); yi=cv2(i); yifit=cv2fit(i);
    %
    %    scatter(log(xi),log(yi))
    %    hold on
    %    scatter(log(xi(1:100)),log(yi(1:100)),'x');
    %    plot(log(xi),log(yifit),'.','markersize',10);
    %    plot(log(xi),log(yifit*chi2inv(0.975,df)./df),'.k');
    %    plot(log(xi),log(yifit*chi2inv(0.025,df)./df),'.k');

    xlabel(hAx1,'Mean expression, log')
    ylabel(hAx1,'CV^2, log')
    if ~isempty(g)
        dt = datacursormode(hFig);        
    end
    hold off

    hAx2=subplot(2,2,2);
    x1 = Xsorted(1,:);
    stem(hAx2, 1:length(x1), x1, 'marker', 'none');
    xlim(hAx2,[1 size(X,2)]);
    title(hAx2, gsorted(1));
    [titxt] = gui.i_getsubtitle(x1);
    subtitle(hAx2, titxt);
    xlabel(hAx2,'Cell Index');
    ylabel(hAx2,'Expression Level');
    
end

dt.UpdateFcn = {@in_myupdatefcn3, g};


    function ChangeAlphaValue(~, ~)
        if h.MarkerFaceAlpha <= 0.05
            h.MarkerFaceAlpha = 1;
        else
            h.MarkerFaceAlpha = h.MarkerFaceAlpha - 0.1;
        end
    end

    function in_HighlightSelectedGenes(~,~)
        %Myc, Oct3/4, Sox2, Klf4
        [glist] = gui.i_selectngenes(SingleCellExperiment(X,g),...
            intersect(upper(g),["MYC", "POU5F1", "SOX2", "KLF4"]));
        if ~isempty(glist)            
            [y,idx]=ismember(glist,g);
            idx=idx(y);            
            % idv = zeros(1, length(hvgidx));
            % idv(idx)=1;
            % h.BrushData = idv;
            for k=1:length(idx)
                dt = datatip(h,'DataIndex',idx(k));
            end
        end
    end

    function HighlightGenes(~, ~)
        idx = zeros(1, length(hvgidx));
        h.BrushData = idx;
        k = gui.i_inputnumk(200, 1, 2000);
        if isempty(k), return; end
        idx(hvgidx(1:k)) = 1;
        h.BrushData = idx;
    end

    function ExportTable(~, ~)
        gui.i_exporttable(T, true, 'Thvgreslist', 'HVGResultTable');
        % Tdegenelist 
        % 'Tviolindata','ViolinPlotTable'
        % 'Thvgreslist', 'HVGResultTable' 
    end

    function ExportGeneNames(~, ~)
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warndlg("No gene is selected.");
            return;
        end
        fprintf('%d genes are selected.\n', sum(ptsSelected));


        gselected=g(ptsSelected);
        [yes,idx]=ismember(gselected,T.genes);
        Tx=T(idx,:);
        Tx=sortrows(Tx,4,'descend');
        if ~all(yes), error('Running time error.'); end
        tgenes=Tx.genes;


        labels = {'Save selected gene names to variable:',...
            'Save HVG table:'};
        vars = {'g','T'};
        values = {tgenes,T};
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
        [yes,idx]=ismember(gselected,T.genes);
        Tx=T(idx,:);
        Tx=sortrows(Tx,4,'descend');
        if ~all(yes), error('Running time error.'); end
        tgenes=Tx.genes;
        
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

    function txt = in_myupdatefcn3(~, event_obj, g)
        if isequal(get(src, 'Parent'), hAx1)
            idx = event_obj.DataIndex;
            txt = {g(idx)};
            x1 = X(idx, :);
            stem(hAx2, 1:length(x1), x1, 'marker', 'none');
            xlim(hAx2,[1 size(X,2)]);
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


% Highly variable genes (HVG) is based on the assumption that genes with
% high variance relative to their mean expression are due to biological
% effects rather than just technical noise. The method seeks to identify
% genes that have a higher variability than expected by considering the
% relationship between variance and mean expression. This relationship is
% difficult to fit, and in practice genes are ranked by their distance
% from a moving median (Kolodziejczyk et al., 2015) or another statistic
% derived from variance is used, e.g. the squared coefficient of variation
% (Brennecke et al. (2013)).
