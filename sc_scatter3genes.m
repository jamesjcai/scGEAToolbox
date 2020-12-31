function sc_scatter3genes(X,genelist,dofit,showdata)
%Scatter3 plot for genes
%

if nargin<2, genelist=[]; end
if nargin<4, showdata=true; end
if nargin<3, dofit=false; end
[lgu,dropr,lgcv,genelist]=sc_genestat(X,genelist);
x=lgu;
y=dropr;
z=lgcv;
if showdata
    scatter3(x,y,z); % 'MarkerEdgeAlpha',.5);
    if ~isempty(genelist)
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn1,genelist};
    end
end
% [xData, yData, zData] = prepareSurfaceData(x,y,z);
% xyz=[xData yData zData]';
if dofit
    [~,xyz1]=sc_splinefit(X,genelist);
%     xyz=[x y z]';
%     % xyz=sortrows([x y z],[1 2])';
%     pieces = 15;
%     s = cumsum([0;sqrt(diff(x(:)).^2 + diff(y(:)).^2 + diff(z(:)).^2)]);
%     pp1 = splinefit(s,xyz,pieces,0.75);
%     xyz1 = ppval(pp1,s);
    hold on
    plot3(xyz1(1,:),xyz1(2,:),xyz1(3,:),'-','linewidth',4);
    % scatter3(xyz1(1,:),xyz1(2,:),xyz1(3,:)); %,'MarkerEdgeAlpha',.8);
end
%grid on
%box on
%legend({'Genes','Spline fit'});
xlabel('Mean, log');
ylabel('Dropout rate (% of zeros)');
zlabel('CV, log');
