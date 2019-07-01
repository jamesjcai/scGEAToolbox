function [X]=norm_bigscale(X)
    [X]=norm_libsize(X);
    X=log2(X+1);
    for k=1:size(X,1)
        if length(unique(X(k,:)))>2
            x_out=shift_values(X(k,:),0,10);
        else
            x_out=zeros(size(X(k,:)));
        end
        X(k,:)=x_out;
    end
	X=X-5;
end


function x_out=shift_values(x,X,Y)
a=min(x);
b=max(x);
x_out=(x-a)/(b-a)*(Y-X)+X;
end

