function sc_explorer2(s1,s2)

% sc_explorer2(s_tsne,s_umap)

explorer2IDX=1:size(s1,1);
assignin('base','explorer2IDX',explorer2IDX);
            
figure;
%https://www.mathworks.com/matlabcentral/answers/153-if-i-have-two-plots-on-the-same-figure-window-how-do-i-use-the-brush-tool-to-highlight-one-data-poi
%https://www.mathworks.com/matlabcentral/answers/385300-how-to-set-the-datasource-of-a-histogram-programmatically

ax1 = subplot(1,2,1);
h1=plot(s1(:,1),s1(:,2),'.');
h1.ZDataSource='explorer2IDX';
ax2 = subplot(1,2,2);
h2=plot(s2(:,1),s2(:,2),'.');
h2.ZDataSource='explorer2IDX';
hLD = linkdata('on');
% hlink = linkprop([ax1,ax2]),{'CameraPosition','CameraUpVector'}); 
% rotate3d on
