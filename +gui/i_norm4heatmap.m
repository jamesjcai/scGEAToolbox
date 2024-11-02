function [Y] = i_norm4heatmap(Y, dim, methodid)

% Y - gene-by-cell matrix
% dim = 2 - by row
if nargin<2, dim = 2; end
if nargin<3, methodid = 1; end

% listitems = {'zscore std', ...
%     'zscore robust', ...
%     'norm 2', ...
%     'norm Inf', ...
%     'scale std', ...
%     'scale mad', ...
%     'scale first', ...
%     'scale iqr', ...
%     'range [0 1]', ...
%     'range [1 10]', ...    
%     'center mean', ...
%     'center median', ...
%     'center scale', ...
%     'medianiqr'};


switch methodid
    case 1
        %Y = zscore(Y, 0, 2);   %         
        Y = normalize(Y, dim, 'zscore');
    case 2
        Y = normalize(Y, dim, 'zscore', 'robust');
        % Y = zscore(Y, 0, 2);
        % Y = reshape(zscore(Y(:)), size(Y));
        % Y = Y ./ (max(abs(Y(:))));
    case 3
        Y = normalize(Y, dim, 'norm', 2);
    case 4
        Y = normalize(Y, dim, 'norm', Inf);
    case 5
        Y = normalize(Y, dim, 'scale', 'std');
    case 6
        Y = normalize(Y, dim, 'scale', 'mad');
    case 7
        Y = normalize(Y, dim, 'scale', 'first');
    case 8
        Y = normalize(Y, dim, 'scale', 'iqr');       
    case 9
        Y = normalize(Y, dim, 'range', [0 1]);
    case 10
        Y = normalize(Y, dim, 'range', [1 10]);
    case 11
        Y = normalize(Y, dim, 'center', 'mean');
    case 12
        Y = normalize(Y, dim, 'center', 'median');
    case 13
        Y = normalize(Y, dim, 'center', 'scale');
    case 14
        Y = normalize(Y, dim, 'medianiqr');
end


        % Y = normalize(Y, dim, "center","median","scale","mad");
        qx = quantile(Y(:), 0.90);
        Y(Y > qx) = qx;
        qx = quantile(Y(:), 0.10);
        %Y(Y<qx)=qx;
        Y(Y < qx) = 0;


