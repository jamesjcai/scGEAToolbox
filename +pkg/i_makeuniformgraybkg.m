

[img,map]=imread('../resources/icon-mat-unfold-less-10 - Copy.gif');
map(255,:)=[0.941176470588235	0.941176470588235	0.941176470588235];
img(img==251)=254;
ptImage= ind2rgb(img, map);
imtool(ptImage)
imwrite(img, map,'../resources/icon-mat-unfold-less-10.gif');





[img,map]=imread('../resources/icon-mat-unfold-more-10 - Copy.gif');
map(255,:)=[0.941176470588235	0.941176470588235	0.941176470588235];
img(img==251)=254;
ptImage= ind2rgb(img, map);
imtool(ptImage)
imwrite(img, map,'../resources/icon-mat-unfold-more-10.gif');


[img,map]=imread('../resources/icon-fa-thumb-tack-10 - Copy.gif');
map(255,:)=[0.941176470588235	0.941176470588235	0.941176470588235];
img(img==251)=254;
ptImage= ind2rgb(img, map);
imtool(ptImage)
imwrite(img, map,'../resources/icon-fa-thumb-tack-10.gif');


[img,map]=imread('../resources/userguiding - Copy.gif');
map(255,:)=[0.941176470588235	0.941176470588235	0.941176470588235];
img(img==251)=254;
ptImage= ind2rgb(img, map);
imtool(ptImage)
imwrite(img, map,'../resources/userguiding.gif');





