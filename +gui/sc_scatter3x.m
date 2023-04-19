function sc_scatter3x(X,Y,genelist,dofit,showdata)
if nargin<3, genelist=[]; end
if nargin<5, showdata=true; end
if nargin<4, dofit=false; end
[lgu1,dropr1,lgcv1,~]=sc_genestat(X,genelist,false);
[lgu2,dropr2,lgcv2,~]=sc_genestat(Y,genelist,false);

if showdata
    scatter3(lgu1,lgcv1,dropr1); % 'MarkerEdgeAlpha',.8);
    hold on
    scatter3(lgu2,lgcv2,dropr2); % 'MarkerEdgeAlpha',.8);
    if ~isempty(genelist)
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcnx,genelist,...
            [lgu1,lgcv1,dropr1],[lgu2,lgcv2,dropr2]};
    end    
end
% [xData, yData, zData] = prepareSurfaceData(x,y,z);
% xyz=[xData yData zData]';
if dofit
    [~,xyz1]=sc_splinefit(X,genelist);
    [~,xyz2]=sc_splinefit(Y,genelist);
    plot3(xyz1(:,1),xyz1(:,2),xyz1(:,3),'-','linewidth',4);
    % scatter3(xyz1(1,:),xyz1(2,:),xyz1(3,:)); %,'MarkerEdgeAlpha',.8);
    plot3(xyz2(:,1),xyz2(:,2),xyz2(:,3),'-','linewidth',4);
end
%grid on
%box on
%legend({'Genes','Spline fit'});
xlabel('Mean, log');
ylabel('CV, log');
zlabel('Dropout rate (% of zeros)');
end