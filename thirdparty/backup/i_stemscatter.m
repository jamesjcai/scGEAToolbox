function i_stemscatter(x,y,z)
if nargin<3
    x=randn(300,1);
    y=randn(300,1);
    z=abs(randn(300,1));
end
if isempty(z)
    warndlg('No expression');
    scatter(x,y,'.');
else
    %scatter(x,y,5);
    stem3(x,y,z,'marker','none','color','m');
    hold on 
    i_myscatter([x y],z);
    hold off
end
%[caz,cel]=view;
view([-45,-45,300]);
end
