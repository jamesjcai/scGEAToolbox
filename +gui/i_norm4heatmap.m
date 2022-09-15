function [Y]=i_norm4heatmap(Y)

methodid=1;
switch methodid
    case 1
        Y=zscore(Y,0,2);
        qx=quantile(Y(:),0.90);
        Y(Y>qx)=qx;
        qx=quantile(Y(:),0.10);
        %Y(Y<qx)=qx;
        Y(Y<qx)=0;
    case 2
        Y=zscore(Y,0,2);
        Y = reshape( zscore(Y(:)), size(Y) );
        Y=Y./(max(abs(Y(:))));
end
