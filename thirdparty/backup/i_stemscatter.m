function i_stemscatter(x,y,z)

stem3(x,y,z,'marker','none');
hold on 
% scatter(x,y,'.');
i_myscatter([x y],z);