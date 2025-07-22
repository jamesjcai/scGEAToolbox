function [x_centre, y_centre,Z_p,a,v,th] = HelixFit(Z,vaxis)
x=Z(1,:);
y=Z(2,:);
z=Z(3,:);
if vaxis=='x'
    [x_centre, y_centre, r, ~] = circleFit ( [y;z] );
    fprintf('vaxis is x default\n');
    a=r;
    [v,th,~,Z_p] = curveFit(x_centre,y_centre,a,x,y,z);
elseif vaxis=='y'
    [x_centre, y_centre, r, ~] = circleFit ( [x;z] );
    fprintf('vaxis is y default\n');
    a=r;
    [v,th,~,Z_p] = curveFit(x_centre,y_centre,a,y,x,z);
elseif vaxis=='z'
    [x_centre, y_centre, r, ~] = circleFit ( [x;y] );
    fprintf('vaxis is z default\n');
    a=r;
    [v,th,~,Z_p] = curveFit(x_centre,y_centre,a,z,x,y);
else
    [x_centre, y_centre, r, sq_error] = circleFit ( [y;z] );
    centre(1,1)=x_centre;
    centre(1,2)=y_centre;
    centre(1,3)=r;
    centre(1,4)=sq_error;
    
    [x_centre, y_centre, r, sq_error] = circleFit ( [x;z] );
    centre(2,1)=x_centre;
    centre(2,2)=y_centre;
    centre(2,3)=r;
    centre(2,4)=sq_error;
    
    [x_centre, y_centre, r, sq_error] = circleFit ( [x;y] );
    centre(3,1)=x_centre;
    centre(3,2)=y_centre;
    centre(3,3)=r;
    centre(3,4)=sq_error;

    min_index=find(centre(:,4)==min(centre(:,4)));
    if min_index==1
       fprintf('vaxis is x\n');
       a=centre(min_index,3);
       [v,th,~,Z_p] = curveFit(x_centre,y_centre,a,x,y,z);
    elseif min_index==2
        fprintf('vaxis is y\n');
        a=centre(min_index,3);
        %[v,th,resnorm,Z_p] = curveFit(a,y,x,z);
        [v,th,~,Z_p] = curveFit(x_centre,y_centre,a,y,x,z);
    elseif min_index==3
        fprintf('vaxis is z\n');
        a=centre(min_index,3);
        %[v,th,resnorm,Z_p] = curveFit(a,z,x,y);
        [v,th,~,Z_p] = curveFit(x_centre,y_centre,a,z,x,y);
    end
end
end
