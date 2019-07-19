function x_out=shift_values(x,X,Y)

a=min(x);
b=max(x);
x_out=(x-a)/(b-a)*(Y-X)+X;



end

