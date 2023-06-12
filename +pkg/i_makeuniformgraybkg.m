

[img,map]=imread('../resources/icon-mat-unfold-less-10 - Copy.gif');
%map(255,:)=1-[0.94 0.94 0.94];
img(img==251)=202;
ptImage= ind2rgb(img, map);
imtool(ptImage)
imwrite(img, map,'../resources/icon-mat-unfold-less-10.gif');


[img,map1]=imread('../resources/icon-mat-unfold-more-10 - Copy.gif');
%map(255,:)=1-[0.94 0.94 0.94];
img(img==251)=202;
ptImage= ind2rgb(img, map);
imtool(ptImage)
imwrite(img, map,'../resources/icon-mat-unfold-more-10.gif');

[img,map1]=imread('../resources/icon-fa-thumb-tack-10 - Copy.gif');
%map(255,:)=1-[0.94 0.94 0.94];
img(img==251)=202;
ptImage= ind2rgb(img, map);
imtool(ptImage)
imwrite(img, map,'../resources/icon-fa-thumb-tack-10.gif');


