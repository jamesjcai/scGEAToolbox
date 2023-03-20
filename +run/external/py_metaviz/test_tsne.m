load fisheriris
rng default % for reproducibility
Y = tsne(meas,'Distance','mahalanobis');
figure;
gscatter(Y(:,1),Y(:,2),species)

D=pdist2(meas,meas,"mahalanobis");
rng default
Y2 = tsne(D,'Distance',@ExampleDistFunc,'NumDimensions',2);
figure;
gscatter(Y2(:,1),Y2(:,2),species);

rng default
[Y3]=pkg.e_tsnebyd(D,2);
figure;
gscatter(Y3(:,1),Y3(:,2),species)
