function hLink=sc_explorer3(s1)

% sc_explorer2(s_tsne,s_umap)

explorer2IDX=s1(:,3);
assignin('base','explorer2IDX',explorer2IDX);


figure;
%https://www.mathworks.com/matlabcentral/answers/153-if-i-have-two-plots-on-the-same-figure-window-how-do-i-use-the-brush-tool-to-highlight-one-data-poi
%https://www.mathworks.com/matlabcentral/answers/385300-how-to-set-the-datasource-of-a-histogram-programmatically

ax1 = subplot(1,2,1);
h1=scatter3(s1(:,1),s1(:,2),s1(:,3),10);
xlabel('dim 1')
ylabel('dim 2')
zlabel('dim 3')
grid on
box on
h1.ZDataSource='explorer2IDX';

ax2 = subplot(1,2,2);
h2=scatter3(s1(:,3),s1(:,4),s1(:,5),10);
xlabel('dim 3')
ylabel('dim 4')
zlabel('dim 5')
grid on
box on
h2.XDataSource='explorer2IDX';
hLD = linkdata('on');


% figure
% ax1 = subplot(2,1,1);
% [X1,Y1,Z1] = peaks;
% surf(X1,Y1,Z1)
% 
% ax2 = subplot(2,1,2);
% [X2,Y2,Z2] = peaks(10);
% surf(X2,Y2,Z2)

%assignin('base','ax1',ax1);
%assignin('base','ax2',ax2);


%hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'}); 
%rotate3d on

%h=findobj(gcf,'type','axes');
%hlink = linkprop(h,{'CameraPosition','CameraUpVector'}); 
%rotate3d on

evalin('base','h=findobj(gcf,''type'',''axes'');');
evalin('base','hlink = linkprop(h,{''CameraPosition'',''CameraUpVector''});');
evalin('base','rotate3d on');

