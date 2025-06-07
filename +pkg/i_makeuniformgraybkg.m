i_aaax('icon-mat-unfold-less-10')
i_aaax('icon-mat-unfold-more-10')
% xxxx
% [img,map]=imread('../assets/icon-mat-unfold-less-10 - Copy.gif');
% map(255,:)=[0.9400 0.9400 0.9400];
% img(img==251)=254;
% ptImage= ind2rgb(img, map);
% imageViewer(ptImage)
% imwrite(img, map,'../assets/icon-mat-unfold-less-10.gif');
%
%
% [img,map]=imread('../assets/icon-mat-unfold-more-10 - Copy.gif');
% map(255,:)=[0.9400 0.9400 0.9400];
% img(img==251)=254;
% ptImage= ind2rgb(img, map);
% imageViewer(ptImage)
% imwrite(img, map,'../assets/icon-mat-unfold-more-10.gif');


[img, map] = imread('../assets/icon-fa-thumb-tack-10 - Copy.gif');
map(255, :) = [0.9400, 0.9400, 0.9400];
img(img == 251) = 254;
ptImage = ind2rgb(img, map);
imageViewer(ptImage)
imwrite(img, map, '../assets/icon-fa-thumb-tack-10.gif');


[img, map] = imread('../assets/userguiding - Copy.gif');
map(255, :) = [0.9400, 0.9400, 0.9400];
img(img == 251) = 254;
ptImage = ind2rgb(img, map);
imageViewer(ptImage)
imwrite(img, map, '../assets/userguiding.gif');


function i_aaax(infilename)
[img, map] = imread(sprintf('../assets/%s - Copy.gif', infilename));
map(255, :) = [0.9400, 0.9400, 0.9400];
img(img == 251) = 254;
ptImage = ind2rgb(img, map);
imageViewer(ptImage)
imwrite(img, map, sprintf('../assets/%s.gif', infilename));
end
