function [v,th,resnorm,Z_p] = curveFit(x_centre,y_centre,a,x,y,z)
v=1;
t=y+z;
f=@(v,x)a*cos(x/v)+x_centre+a*sin(x/v)+y_centre;
%使用lsqcurvefit
[A,resnorm]=lsqcurvefit(f,v,x,t);
v=A;
th=x./v;
y_new=a*cos(th)+x_centre;
z_new=a*sin(th)+y_centre;
Z_p=[x;y_new;z_new];
end

