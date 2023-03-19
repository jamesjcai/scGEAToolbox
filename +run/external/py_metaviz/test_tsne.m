load fisheriris
rng default % for reproducibility
Y = tsne(meas,'Distance','euclidean');
%gscatter(Y(:,1),Y(:,2),species)

D=pdist2(meas,meas);
Y2 = tsne(D,'Distance',@ExampleDistFunc);
gscatter(Y2(:,1),Y2(:,2),species)

