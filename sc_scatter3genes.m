function sc_scatter3genes(X,g,dofit,showdata)
%Scatter3 plot for genes

if nargin<4, showdata=true; end
if nargin<3, dofit=true; end
if nargin<2, g=[]; end
[lgu,dropr,lgcv,g]=sc_genestat(X,g);
x=lgu;
y=dropr;
z=lgcv;
if showdata
    FigureHandle=figure;
    hAx = axes('Parent', FigureHandle);
    UitoolbarHandle = uitoolbar('Parent', FigureHandle);
    set(UitoolbarHandle, 'Tag', 'FigureToolBar', ...
        'HandleVisibility', 'off', 'Visible', 'on');

    ptlabelclusters = uipushtool(UitoolbarHandle);
    [img, map] = imread(fullfile(matlabroot, ...
        'toolbox', 'matlab', 'icons', 'plotpicker-scatter.gif'));
    ptImage = ind2rgb(img, map);
    ptlabelclusters.CData = ptImage;
    ptlabelclusters.Tooltip = 'Label clusters';
    ptlabelclusters.ClickedCallback = @ExportGeneNames;


    ptlabelclusters = uipushtool(UitoolbarHandle);
    [img, map] = imread(fullfile(matlabroot, ...
        'toolbox', 'matlab', 'icons', 'plotpicker-scatter.gif'));
    ptImage = ind2rgb(img, map);
    ptlabelclusters.CData = ptImage;
    ptlabelclusters.Tooltip = 'Label clusters';
    ptlabelclusters.ClickedCallback = @HighlightGenes;


    h=scatter3(hAx,x,y,z);  % 'filled','MarkerFaceAlpha',.5);
    if ~isempty(g)
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn1,g};
    end
end
% [xData, yData, zData] = prepareSurfaceData(x,y,z);
% xyz=[xData yData zData]';
if dofit
    [~,xyz1]=sc_splinefit(X,g);
%     xyz=[x y z]';
%     % xyz=sortrows([x y z],[1 2])';
%     pieces = 15;
%     s = cumsum([0;sqrt(diff(x(:)).^2 + diff(y(:)).^2 + diff(z(:)).^2)]);
%     pp1 = splinefit(s,xyz,pieces,0.75);
%     xyz1 = ppval(pp1,s);
    hold on
    plot3(xyz1(1,:),xyz1(2,:),xyz1(3,:),'-','linewidth',4);
    % scatter3(xyz1(1,:),xyz1(2,:),xyz1(3,:)); %,'MarkerEdgeAlpha',.8);

    [~,d]=dsearchn(xyz1.',[x y z]);
    [~,hvgidx]=sort(d,'descend');
    
    %g(idx20)
    
end

%grid on
%box on
%legend({'Genes','Spline fit'});
xlabel('Mean, log');
ylabel('Dropout rate (% of zeros)');
zlabel('CV, log');

   function HighlightGenes(~,~)
        %h.MarkerIndices=idx20;
        k=gui.i_inputnumk(50);
        if isempty(k), return; end
        idx=zeros(1,length(hvgidx));
        idx(hvgidx(1:k))=1;
        h.BrushData=idx;
        % datatip(h, 'DataIndex', idx20);
        %h2=scatter3(x(idx20),y(idx20),z(idx20),'rx');  % 'filled','MarkerFaceAlpha',.5);
    end

function ExportGeneNames(~,~)
        
        ptsSelected = logical(h.BrushData.');
        if ~any(ptsSelected)
            warndlg("No gene is selected.");
            return;
        end
        fprintf('%d genes are selected.\n',sum(ptsSelected));

        labels = {'Save gene names to variable:'}; 
        vars = {'g'};
        values = {g(ptsSelected)};
        export2wsdlg(labels,vars,values,...
                     'Save Data to Workspace');
end
end