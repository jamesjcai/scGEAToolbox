tissue_paths='spatial/tissue_positions_list.csv';
image_paths='spatial/tissue_lowres_image.png';
img=imread(image_paths);
% imtool(a);
T=readtable(tissue_paths);
assert(all(ismember(sce.c_cell_id,string(T.Var1))))
[~,idx]=ismember(sce.c_cell_id,string(T.Var1));
t=T(idx,:);

%%
s=[t.Var3,t.Var4];
scgeatool(sce.X,sce.g,s,sce.c)
hold on
xImage = [0 60; 0 60];   % The x data for the image corners
zImage = [-1 -1; -1 -1];   % The z data for the image corners
yImage = [140 140; 0 0];   % The z data for the image corners

img=flip(img);
img=flip(img,2);
%img=imrotate(img,-180);
surf(xImage,yImage,zImage,...    % Plot the surface
     'CData',img,...
     'FaceColor','texturemap');
view(2);
