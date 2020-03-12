function sc_stemscatter(x,y,z)
if nargin<3
    x=randn(300,1);
    y=randn(300,1);
    z=abs(randn(300,1));
end
stem3(x,y,z,'marker','none');
hold on 
% scatter(x,y,'.');
i_myscatter([x y],z);
hold off
