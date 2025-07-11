function [xyz, xyz1, gsorted] = sc_scatter3genes(X, g, dofit, showdata, ...
    parentfig)
%Scatter3 plot for genes

if nargin < 5, parentfig = []; end
if nargin < 4, showdata = true; end
if nargin < 3, dofit = true; end
if nargin < 2 || isempty(g), g = string(1:size(X,1)); end
xyz=[]; xyz1=[];

[lgu, dropr, lgcv, gsorted, Xsorted] = sc_genestat(X, g);

x = lgu;
y = lgcv;
z = dropr;

fw = gui.myWaitbar(parentfig);


hx=gui.myFigure(parentfig);
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
    %h=scatter3(hAx,x,y,z);  % 'filled','MarkerFaceAlpha',.5);
    hAx1 = subplot(2,2,[1 3]);
    h = scatter3(hAx1, x, y, z, 'filled', 'MarkerFaceAlpha', .1);

    if ~isempty(gsorted)
        dt = datacursormode(hFig);
        % dt.UpdateFcn = {@i_myupdatefcn1, g};
    else
        dt = [];
    end
end    
    %grid on
    %box on
    %legend({'Genes','Spline fit'});
    xlabel(hAx1,'Mean+1, log');
    ylabel(hAx1,'CV+1, log');
    zlabel(hAx1,'Dropout rate (% of zeros)');

        
% [xData, yData, zData] = prepareSurfaceData(x,y,z);
% xyz=[xData yData zData]';

if dofit
    try
        [~, ~, ~, xyz1] = sc_splinefit(Xsorted, gsorted);
    catch ME
        rethrow(ME);
    end

    %     xyz=[x y z]';
    %     % xyz=sortrows([x y z],[1 2])';
    %     pieces = 15;
    %     s = cumsum([0;sqrt(diff(x(:)).^2 + diff(y(:)).^2 + diff(z(:)).^2)]);
    %     pp1 = splinefit(s,xyz,pieces,0.75);
    %     xyz1 = ppval(pp1,s);
    hold(hAx1,'on');
    xyz = [x y z];
    %assignin("base","xyz1a", xyz1);
    %assignin("base","xyz1a", xyz);
        

    plot3(hAx1, xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), '-', 'linewidth', 4);
    % scatter3(xyz1(:,1),xyz1(:,2),xyz1(:,3)); %,'MarkerEdgeAlpha',.8);

    [~, d] = dsearchn(xyz1, [x, y, z]);

    fitmeanv=xyz1(:,1);
    d(x>max(fitmeanv))=d(x>max(fitmeanv))./100;
    d(x<min(fitmeanv))=d(x<min(fitmeanv))./10;
    d((y-xyz1(:, 2))<0)=d((y-xyz1(:, 2))<0)./100;

    [sortedd, hvgidx] = sort(d, 'descend');

    hvg=gsorted(hvgidx);
    lgu=lgu(hvgidx);
    lgcv=lgcv(hvgidx);
    dropr=dropr(hvgidx);    
    
    gene=hvg(:);
    T=table(gene,sortedd,hvgidx,lgu,lgcv,dropr);
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


hAx2 = subplot(2,2,2);
x1=Xsorted(hvgidx(1),:);
%sh = stem(hAx2, 1:length(x1), x1, 'marker', 'none');

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


    % function in_ShowProfile(~, ~)
    %     idx = 1;
    %     x1 = X(idx, :);
    %     stem(hAx2, 1:length(x1), x1, 'marker', 'none');
    %     xlim(hAx2,[1 size(X,2)]);
    %     title(hAx2, g(idx));
    %     xlabel(hAx2,'Cell Index')
    %     ylabel(hAx2,'Expression Level')
    % end

    function ChangeAlphaValue(~, ~)
        if h.MarkerFaceAlpha <= 0.05
            h.MarkerFaceAlpha = 1;
        else
            h.MarkerFaceAlpha = h.MarkerFaceAlpha - 0.1;
        end
    end            

    function in_callback_HighlightTopHVGs(~, ~)
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

    function in_callback_ExportTable(~, ~)                
        gui.i_exporttable(T, true, 'Tsplinefitg', 'SplinefitGTable');
    end

    function in_HighlightSelectedGenes(~,~,typeid)        
        if nargin<3, typeid = 1; end

           switch typeid
               case 1
                    %Myc, Oct3/4, Sox2, Klf4
                    [glist] = gui.i_selectngenes(SingleCellExperiment(Xsorted,gsorted),...
                        intersect(upper(gsorted),["MYC", "POU5F1", "SOX2", "KLF4"]));                    
               case 2
                    gsorted = T.(T.Properties.VariableNames{1});
                    
                   if gui.i_isuifig(parentfig)
                        [indx2, tf2] = gui.myListdlg(parentfig, gsorted, 'Select genes:');
                    else
                        [indx2, tf2] = listdlg('PromptString', ...
                            'Select genes:', ...
                            'SelectionMode', 'multiple', ...
                            'ListString', gsorted, ...
                            'ListSize', [220, 300]);
                    end


                    if tf2 == 1
                        glist = gsorted(indx2);
                    else
                        return;
                    end
           end

        if ~isempty(glist)            
            [yes,idx]=ismember(glist,gsorted);
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
            gui.myWarndlg(parentfig, "No gene is selected.");
            return;
        end
        fprintf('%d genes are selected.\n', sum(ptsSelected));

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
            gui.myWarndlg(parentfig, "No gene is selected.");
            return;
        end
        fprintf('%d genes are selected.\n', sum(ptsSelected));        

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

            % figure; stem(1:length(x1), x1, 'marker', 'none');

            if ~isempty(sh) && isvalid(sh)
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
