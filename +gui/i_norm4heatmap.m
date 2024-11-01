function [Y] = i_norm4heatmap(Y, dim)

% Y - gene-by-cell matrix
% dim = 2 - by row
if nargin<2, dim = 2; end

methodid = 1;
switch methodid
    case 1
        %Y = zscore(Y, 0, 2);   % 
        
        Y = normalize(Y, dim, 'zscore');

        % Y = normalize(Y, dim, "center","median","scale","mad");
        qx = quantile(Y(:), 0.90);
        Y(Y > qx) = qx;
        qx = quantile(Y(:), 0.10);
        %Y(Y<qx)=qx;
        Y(Y < qx) = 0;
    case 2
        Y = zscore(Y, 0, 2);
        Y = reshape(zscore(Y(:)), size(Y));
        Y = Y ./ (max(abs(Y(:))));
end
