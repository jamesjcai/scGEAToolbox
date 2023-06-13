i_aaax('icon-mat-unfold-less-10')
i_aaax('icon-mat-unfold-more-10')

% [img,map]=imread('../resources/icon-mat-unfold-less-10 - Copy.gif');
% map(255,:)=[0.9400 0.9400 0.9400];
% img(img==251)=254;
% ptImage= ind2rgb(img, map);
% imtool(ptImage)
% imwrite(img, map,'../resources/icon-mat-unfold-less-10.gif');
% 
% 
% [img,map]=imread('../resources/icon-mat-unfold-more-10 - Copy.gif');
% map(255,:)=[0.9400 0.9400 0.9400];
% img(img==251)=254;
% ptImage= ind2rgb(img, map);
% imtool(ptImage)
% imwrite(img, map,'../resources/icon-mat-unfold-more-10.gif');


[img,map]=imread('../resources/icon-fa-thumb-tack-10 - Copy.gif');
map(255,:)=[0.9400 0.9400 0.9400];
img(img==251)=254;
ptImage= ind2rgb(img, map);
imtool(ptImage)
imwrite(img, map,'../resources/icon-fa-thumb-tack-10.gif');


[img,map]=imread('../resources/userguiding - Copy.gif');
map(255,:)=[0.9400 0.9400 0.9400];
img(img==251)=254;
ptImage= ind2rgb(img, map);
imtool(ptImage)
imwrite(img, map,'../resources/userguiding.gif');



function i_aaax(infilename)
    [img,map]=imread(sprintf('../resources/%s - Copy.gif',infilename));
    map(255,:)=[0.9400 0.9400 0.9400];
    img(img==251)=254;
    ptImage= ind2rgb(img, map);
    imtool(ptImage)
    imwrite(img, map,sprintf('../resources/%s.gif',infilename));
end

